import json
import logging
import os
import shutil
import sys
import time
from io import StringIO
from pathlib import Path
from typing import Mapping, Optional, Union

from Bio import Align
from Bio.PDB import MMCIFIO, MMCIFParser
from colorama import Fore, Style

logger = logging.getLogger("logger")


# 색상 로깅을 위한 사용자 정의 포맷터
class ColoredFormatter(logging.Formatter):
    # 각 로그 수준에 대한 색상 코드 정의
    LEVEL_COLORS = {
        logging.DEBUG: Fore.BLUE,
        logging.INFO: Fore.WHITE,
        logging.WARNING: Fore.YELLOW,
        logging.ERROR: Fore.RED,
        logging.CRITICAL: Fore.RED + Style.BRIGHT,
    }

    def format(self, record):
        # 로그 수준에 대한 색상 가져오기
        level_color = self.LEVEL_COLORS.get(record.levelno, "")
        # 로그 메시지 포맷 지정
        formatted_message = super().format(record)
        # 색상이 추가된 메시지 반환
        return f"{level_color}{formatted_message}{Style.RESET_ALL}"


# 로깅 설정
def setup_logger():
    logger = logging.getLogger("logger")
    logger.setLevel(logging.DEBUG)  # 최소 로깅 수준 설정

    # 스트림 핸들러 생성 (콘솔 출력)
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)

    # 사용자 정의 포맷터 설정
    formatter = ColoredFormatter("%(asctime)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)

    # 로거에 핸들러 추가
    logger.addHandler(handler)

    return logger


def get_chains(mmcif_file):
    """MMCIF 파일의 체인 목록을 반환합니다."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("template", mmcif_file)
    chains = []
    for model in structure:
        for chain in model:
            chains.append(chain.id)
    return chains


def extract_sequence_from_mmcif(mmcif_file):
    """MMCIF 파일에서 서열을 추출합니다."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("template", mmcif_file)
    sequence = ""
    model = structure[0]  # 단일 모델/체인만 있다고 가정
    for chain in model:
        for residue in chain:
            if residue.id[0] == " ":  # 헤테로 원자 제외
                sequence += residue.resname[0]  # 첫 글자만 가져오도록 단순화
    return sequence


# 코드는 https://github.com/google-deepmind/alphafold3 에서 가져옴
def query_to_hit_mapping(
    query_aligned: str, template_aligned: str
) -> Mapping[int, int]:
    """0-기반 쿼리 인덱스를 히트 인덱스로 매핑합니다."""
    query_to_hit_mapping_out = {}
    hit_index = 0
    query_index = 0
    for q_char, t_char in zip(query_aligned, template_aligned):
        # 템플릿에 갭 삽입됨
        if q_char == "-":
            query_index += 1
        # 템플릿에서 삭제된 잔기 (쿼리에서는 갭).
        elif t_char == "-":
            hit_index += 1
        # 쿼리와 템플릿 모두에 있는 정상적으로 정렬된 잔기. 매핑에 추가.
        else:
            query_to_hit_mapping_out[query_index] = hit_index
            query_index += 1
            hit_index += 1
    return query_to_hit_mapping_out


def align_and_map(query_seq, template_seq):
    """두 서열을 정렬하고 인덱스를 매핑합니다."""
    # 쌍별 정렬 수행
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(query_seq, template_seq)
    alignment = alignments[0]  # 최적의 정렬 선택

    formatted_alignment = alignment._format_generalized().replace(" ", "")
    query_aligned, _, template_aligned, _ = formatted_alignment.split("\n")

    # 정렬된 서열 매핑
    aligned_mapping = query_to_hit_mapping(query_aligned, template_aligned)

    query_indices = []
    template_indices = []
    for template_index, query_index in aligned_mapping.items():
        query_indices.append(query_index)
        template_indices.append(template_index)

    return query_indices, template_indices


def get_mmcif(
    cif,
    pdb_id,
    chain_id,
    start,
    end,
    tmpdir=None,
):
    """
    CIF 파일에서 체인을 추출하고 지정된 체인, 잔기 및 메타데이터만 포함하는
    새로운 CIF 문자열을 반환합니다.
    """

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, cif)

    # CIF 파일에서 릴리스 날짜 추출
    mmcif_dict = parser._mmcif_dict
    headers_to_keep = [
        "_entry.id",
        "_entry.title",
        "_entry.deposition_date",
        "_pdbx_audit_revision_history.revision_date",
    ]
    filtered_metadata = {
        key: mmcif_dict[key] for key in headers_to_keep if key in mmcif_dict
    }

    # 누락된 경우 메타데이터 생성
    if "_pdbx_audit_revision_history.revision_date" not in filtered_metadata:
        filtered_metadata["_pdbx_audit_revision_history.revision_date"] = time.strftime(
            "%Y-%m-%d"
        )

    # 다중 모델 템플릿(예: NMR)의 경우 단일 대표 모델 선택
    # 구조에 첫 번째 모델의 복사본 추가

    if len(structure) > 1:
        for model_index in range(1, len(structure)):
            structure.detach_child(structure[model_index].get_id())

    for model in structure:
        chain_to_del = []
        for chain in model:
            if chain.id != chain_id:
                chain_to_del.append(chain.id)
                continue

        for unwanted_chain in chain_to_del:
            model.detach_child(unwanted_chain)

        for chain in model:
            res_to_del = []
            for i, res in enumerate(chain):
                rel_pos = i + 1
                if rel_pos < start or rel_pos > end or res.id[0] != " ":
                    res_to_del.append(res)

            for res in res_to_del:
                chain.detach_child(res.id)

    # 필터링된 구조를 새 CIF 파일에 저장
    io = MMCIFIO()
    io.set_structure(structure)
    filtered_output = (
        f"{pdb_id}_{chain_id}.cif"
        if tmpdir is None
        else f"{tmpdir}/{pdb_id}_{chain_id}.cif"
    )
    io.save(filtered_output)

    # 필터링된 구조를 파싱하여 메타데이터가 없는 수정된 MMCIF 가져오기
    structure = parser.get_structure(pdb_id, filtered_output)
    mmcif_dict = parser._mmcif_dict

    # 필터링된 메타데이터를 MMCIF 딕셔너리에 추가
    mmcif_dict.update(filtered_metadata)

    # 원하는 메타데이터가 포함된 수정된 MMCIF를 문자열로 저장
    string_io = StringIO()
    io.set_dict(mmcif_dict)
    io.save(string_io)

    os.unlink(filtered_output)

    return string_io.getvalue()


# 입력 JSON의 각 서열에 대해 실행
def get_custom_template(
    sequence,
    target_id,
    custom_template,
    custom_template_chain,
):
    """입력 JSON에 사용자 정의 템플릿을 추가합니다."""

    # 여기에 사용자 정의 템플릿을 실행하는 코드를 추가하되, 결과는 출력으로 옮김
    # 기존 템플릿이 있는 경우 유지
    if "templates" not in sequence["protein"]:
        templates = []
    else:
        templates = sequence["protein"]["templates"]

    input_sequence = sequence["protein"]["sequence"]
    seq_id = sequence["protein"]["id"]
    if isinstance(seq_id, list):
        if target_id and target_id not in seq_id:
            return sequence
    if isinstance(seq_id, str):
        if target_id and target_id != seq_id:
            return sequence

    if not os.path.exists(custom_template):
        msg = f"사용자 정의 템플릿 파일 {custom_template}을(를) 찾을 수 없습니다."
        logger.critical(msg)
        raise FileNotFoundError()

    chain_info = get_chains(custom_template)
    if len(chain_info) != 1 and not custom_template_chain:
        msg = f"사용자 정의 템플릿 파일 {custom_template}에 \
{len(chain_info)}개의 체인이 포함되어 있습니다. --custom_template_chain으로 사용할 체인을 지정하십시오."
        raise ValueError(msg)

    if custom_template_chain and custom_template_chain not in chain_info:
        msg = f"사용자 정의 템플릿 파일 {custom_template}에 \
체인 {custom_template_chain}이(가) 포함되어 있지 않습니다."
        raise ValueError(msg)

    if not custom_template_chain:
        custom_template_chain = chain_info[0]

    template = {}
    cif_str = get_mmcif(
        custom_template,
        "custom",
        custom_template_chain,
        1,
        len(input_sequence),
    )

    template["mmcif"] = cif_str
    template_seq = extract_sequence_from_mmcif(StringIO(cif_str))
    query_indices, template_indices = align_and_map(input_sequence, template_seq)

    template["queryIndices"] = query_indices
    template["templateIndices"] = template_indices

    # 사용자 정의 템플릿을 템플릿 목록의 시작 부분에 추가
    templates.insert(0, template)

    # JSON에 템플릿 추가
    sequence["protein"]["templates"] = templates

    # 출력 JSON 저장
    return sequence


def make_dir(dir_path: Union[str, Path], overwrite: bool = False):
    """
    디렉터리를 만들고 Path 객체를 반환합니다.

    Args:
        dir_path: 만들 디렉터리 경로입니다.
        overwrite: 디렉터리가 이미 있는 경우 삭제할지 여부입니다.

    Returns:
        생성된 디렉터리에 대한 Path 객체입니다.
    """
    dir_path = Path(dir_path)
    if dir_path.exists():
        if overwrite:
            shutil.rmtree(dir_path)
        else:
            logger.error(
                f"디렉터리 {dir_path}이(가) 이미 존재합니다. 바꾸려면 --override를 사용하십시오."
            )
            raise FileExistsError()

    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path


def check_input_json(
    input_json: Union[str, Path],
    output_dir: Optional[Union[str, Path]] = None,
    use_af3_templates: bool = False,
    test: bool = False,
):
    """
    입력 JSON 파일에서 누락된 필드를 확인하고 기본값을 추가합니다.

    Args:
        input_json: 입력 JSON 파일 경로입니다.
        output_dir: 출력 디렉터리 경로입니다.
        use_af3_templates: AlphaFold3 템플릿을 사용할지 여부입니다.
        test: 테스트 모드인지 여부입니다.

    Returns:
        Path: 처리된 JSON 파일의 경로입니다.
    """
    input_json = Path(input_json)

    output_json = (
        input_json.parent.joinpath("abc_" + input_json.name)
        if output_dir is None
        else Path(output_dir).joinpath("abc_" + input_json.name)
    )
    with open(input_json, "r") as f:
        input_data = json.load(f)

    for sequence in input_data["sequences"]:
        for sequence_type in sequence:
            if "unpairedMsaPath" in sequence[sequence_type]:
                msa_path = sequence[sequence_type]["unpairedMsaPath"]
                if not os.path.exists(msa_path):
                    logger.error(f"MSA 파일 {msa_path}을(를) 찾을 수 없습니다.")
                    sys.exit(1)
                with open(msa_path, "r") as f:
                    msa = f.read()
                sequence[sequence_type]["unpairedMsa"] = msa
                del sequence[sequence_type]["unpairedMsaPath"]
            if "unpairedMsa" in sequence[sequence_type]:
                if "templates" not in sequence[sequence_type]:

                    sequence[sequence_type]["templates"] = (
                        None if use_af3_templates else []
                    )

                if "pairedMsa" not in sequence[sequence_type]:
                    sequence[sequence_type]["pairedMsa"] = ""

    if test:
        return input_data
    with open(output_json, "w") as f:
        json.dump(input_data, f, indent=4)

    return output_json


def make_dummy_af3_db(output_dir):
    """AlphaFold3용 더미 데이터베이스를 생성합니다."""
    dummy_af3_db = output_dir.joinpath("af3_db")
    dummy_files = [
        "bfd-first_non_consensus_sequences.fasta",
        "rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta",
        "uniprot_all_2021_04.fa",
        "mgy_clusters_2022_05.fa",
        "nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta",
        "pdb_seqres_2022_09_28.fasta",
        "rnacentral_active_seq_id_90_cov_80_linclust.fasta",
        "uniref90_2022_05.fa"
        ]
    dummy_dirs = ["mmcif_files"]

    dummy_af3_db.mkdir(parents=True, exist_ok=True)
    for dummy_file in dummy_files:
        with open(dummy_af3_db.joinpath(dummy_file), "w") as f:
            f.write("")
    for dummy_dir in dummy_dirs:
        dummy_af3_db.joinpath(dummy_dir).mkdir(parents=True, exist_ok=True)

    return dummy_af3_db
