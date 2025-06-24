#!/usr/bin/env python

import gzip
import json
import logging
import math
import os
import random
import shutil
import subprocess
import tarfile
import tempfile
import time
from io import StringIO
from pathlib import Path
from typing import List, Sequence, Union

import pandas as pd
import requests  # type: ignore
from tqdm.autonotebook import tqdm

from abcfold.argparse_utils import (custom_template_argpase_util,
                                    mmseqs2_argparse_util)
from abcfold.scripts.abc_script_utils import (align_and_map,
                                              extract_sequence_from_mmcif,
                                              get_custom_template, get_mmcif)

logger = logging.getLogger("logger")

TQDM_BAR_FORMAT = (
    "{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]"
)

MODULE_OUTPUT_POS = {
    "align":        4,
    "convertalis":  4,
    "expandaln":    5,
    "filterresult": 4,
    "lndb":         2,
    "mergedbs":     2,
    "mvdb":         2,
    "pairaln":      4,
    "result2msa":   4,
    "search":       3,
}


class MMseqs2Exception(Exception):
    """MMseqs2 API 관련 오류를 위한 사용자 정의 예외 클래스입니다."""
    def __init__(self):

        msg = "MMseqs2 API에서 오류가 발생하고 있습니다. 입력이 유효한 단백질 서열인지 확인하십시오. \
오류가 지속되면 1시간 후에 다시 시도하십시오."
        logger.error(msg)
        super().__init__(msg) # Pass the message to the base Exception class


def add_msa_to_json(
    input_json,
    mmseqs_db,
    templates,
    num_templates,
    chai_template_output,
    custom_template,
    custom_template_chain,
    target_id,
    input_params=None,
    output_json=None,
    to_file=True,
):
    """
    AlphaFold3 입력 JSON 객체에 MMseqs2 MSA 및/또는 사용자 정의 템플릿을 추가합니다.

    Args:
        input_json (str or Path): 입력 AlphaFold3 JSON 파일 경로.
        mmseqs_db (str or Path or None): 로컬 MMseqs2 데이터베이스 경로.
                                         None이면 MMseqs2 웹 서버를 사용합니다.
        templates (bool): 템플릿을 검색할지 여부.
        num_templates (int): 사용할 템플릿 수.
        chai_template_output (str or Path or bool): Chai 템플릿 히트를 저장할 경로,
                                                   또는 사용하지 않는 경우 False.
        custom_template (list or None): 사용자 정의 템플릿 CIF 파일 경로 목록.
        custom_template_chain (list or None): 사용자 정의 템플릿에 대한 체인 ID 목록.
        target_id (list or None): 사용자 정의 템플릿에 대한 대상 ID 목록.
        input_params (dict, optional): 미리 로드된 입력 JSON(dict). 기본값은 None입니다.
        output_json (str or Path, optional): 수정된 JSON을 저장할 경로.
                                            None이고 to_file이 True이면
                                            "_mmseqs" 접미사가 붙은 새 파일을 만듭니다.
        to_file (bool, optional): 출력을 파일에 저장할지 여부. 기본값은 True입니다.

    Returns:
        dict: 수정된 AlphaFold3 JSON 객체.

    Raises:
        FileNotFoundError: 사용자 정의 템플릿 파일을 찾을 수 없는 경우.
        ValueError: 템플릿 인수에 불일치가 있는 경우.
        MMseqs2Exception: MMseqs2 API가 오류를 반환하는 경우.
    """
    if input_params is None:
        with open(input_json, "r") as f:
            input_params = json.load(f)

    for sequence in input_params["sequences"]:
        if "protein" in sequence:
            input_id = sequence["protein"]["id"]
            input_sequence = sequence["protein"]["sequence"]
            with tempfile.TemporaryDirectory() as tmpdir:
                if mmseqs_db:
                    logger.info(f"로컬 MMseqs2를 서열에 실행 중: {input_sequence}")
                    if templates:
                        a3m_lines, templates = run_local_mmseqs(
                            input_sequence,
                            Path(tmpdir),
                            use_templates=True,
                            num_templates=num_templates,
                            mmseqs_db=Path(mmseqs_db),
                        )
                    else:
                        a3m_lines = run_local_mmseqs(
                            input_sequence,
                            Path(tmpdir),
                            use_templates=False,
                            mmseqs_db=Path(mmseqs_db),
                        )
                else:
                    logger.info(f"MMseqs2를 서열에 실행 중: {input_sequence}")
                    # MMseqs2를 실행하여 unpaired MSA 가져오기
                    if templates:
                        a3m_lines, templates = run_mmseqs(
                            input_sequence,
                            tmpdir,
                            use_templates=True,
                            num_templates=num_templates,
                        )

                        for i in input_id:
                            table = pd.read_csv(
                                f"{tmpdir}/pdb70.m8",
                                delimiter="\t",
                                header=None,
                                names=[
                                    "query_id",
                                    "subject_id",
                                    "pident",
                                    "length",
                                    "mismatch",
                                    "gapopen",
                                    "query_start",
                                    "query_end",
                                    "subject_start",
                                    "subject_end",
                                    "evalue",
                                    "bitscore",
                                    "comment",
                                ],
                            )

                            table["query_id"] = i

                            if chai_template_output:
                                if os.path.exists(chai_template_output):
                                    table.to_csv(
                                        chai_template_output,
                                        sep="\t",
                                        index=False,
                                        header=False,
                                        mode="a",
                                    )
                                else:
                                    table.to_csv(
                                        chai_template_output,
                                        sep="\t",
                                        index=False,
                                        header=False,
                                    )

                    else:
                        a3m_lines = run_mmseqs(
                            input_sequence,
                            tmpdir,
                            use_templates=False)
                        templates = []

                if custom_template:
                    for template_file in custom_template: # 변수명 변경: template -> template_file (상위 스코프 templates와 충돌 방지)
                        if not os.path.exists(template_file):
                            msg = f"사용자 정의 템플릿 파일 {template_file}을(를) 찾을 수 없습니다."
                            logger.critical(msg)
                            raise FileNotFoundError()
                        # 단백질 서열에만 템플릿을 추가할 수 있으므로 입력 JSON에
                        # 여러 단백질 서열이 있는지 확인합니다.
                        if (
                            len(
                                [
                                    x
                                    for x in input_params["sequences"]
                                    if "protein" in x.keys()
                                ]
                            )
                            > 1
                            and not target_id
                        ):
                            msg = "입력 JSON에서 여러 서열이 발견되었습니다. \
사용자 정의 템플릿을 올바른 서열에 추가하려면 대상 ID를 지정하십시오."
                            raise ValueError(msg)

                    if target_id and len(target_id) > 1:
                        if (len(custom_template) != len(target_id)) or (
                            len(custom_template_chain) != len(target_id)
                        ):
                            msg = "여러 대상에 대한 템플릿을 제공하는 경우, 대상 ID의 수는 \
사용자 정의 템플릿 및 사용자 정의 템플릿 체인의 수와 일치해야 합니다."
                            raise ValueError(msg)
                        custom_templates_zip = zip( # 변수명 변경: custom_templates -> custom_templates_zip
                            target_id, custom_template, custom_template_chain
                        )
                    else:
                        if len(custom_template) != len(custom_template_chain):
                            msg = "사용자 정의 템플릿의 수는 사용자 정의 템플릿 체인의 수와 일치해야 합니다."
                            raise ValueError(msg)
                        # 단일 대상 ID가 제공되면 모든 사용자 정의 템플릿이
                        # 동일한 대상을 위한 것이라고 가정합니다.
                        if target_id:
                            target_ids = [target_id[0]] * len(custom_template)
                        else:
                            target_ids = [None] * len(custom_template)
                        custom_templates_zip = zip( # 변수명 변경: custom_templates -> custom_templates_zip
                            target_ids, custom_template, custom_template_chain
                        )

                    for i in custom_templates_zip: # 변수명 변경: custom_templates -> custom_templates_zip
                        tid, c_tem, c_tem_chn = i
                        sequence = get_custom_template(
                            sequence,
                            tid,
                            c_tem,
                            c_tem_chn,
                        )

                # JSON에 unpaired MSA 추가
                sequence["protein"]["unpairedMsa"] = a3m_lines[0]
                sequence["protein"]["pairedMsa"] = ""
                sequence["protein"]["templates"] = templates

    if to_file:
        if output_json:
            with open(output_json, "w") as f:
                json.dump(input_params, f)
        else:
            # Path 객체로 변환 후 suffix 변경
            input_json_path = Path(input_json)
            output_json = input_json_path.with_name(f"{input_json_path.stem}_mmseqs.json")
            with open(output_json, "w") as f:
                json.dump(input_params, f)

    return input_params


# ColabFold에서 약간 수정된 코드: https://github.com/sokrypton/ColabFold
def run_mmseqs(
    x,
    prefix,
    use_env=True,
    use_filter=True,
    use_templates=False,
    filter=None,
    use_pairing=False,
    host_url="https://a3m.mmseqs.com",
    num_templates=20
) -> Sequence[object]:
    submission_endpoint = "ticket/pair" if use_pairing else "ticket/msa"

    def submit(seqs, mode, N=101):
        n, query = N, ""
        for seq in seqs:
            query += f">{n}\n{seq}\n"
            n += 1

        res = requests.post(
            f"{host_url}/{submission_endpoint}", data={"q": query, "mode": mode}
        )
        try:
            out = res.json()
        except ValueError:
            logger.error(f"서버가 JSON으로 응답하지 않았습니다: {res.text}")
            out = {"status": "ERROR"}
        return out

    def status(ID):
        res = requests.get(f"{host_url}/ticket/{ID}")
        try:
            out = res.json()
        except ValueError:
            logger.error(f"서버가 JSON으로 응답하지 않았습니다: {res.text}")
            out = {"status": "ERROR"}
        return out

    def download(ID, path):
        res = requests.get(f"{host_url}/result/download/{ID}")
        with open(path, "wb") as out_file: # 변수명 out -> out_file (외부 out과 충돌 방지)
            out_file.write(res.content)

    # 입력 x 처리
    seqs = [x] if isinstance(x, str) else x

    # 이전 옵션과의 호환성
    if filter is not None:
        use_filter = filter

    # 모드 설정
    if use_filter:
        mode = "env" if use_env else "all"
    else:
        mode = "env-nofilter" if use_env else "nofilter"

    if use_pairing:
        mode = ""
        use_templates = False
        use_env = False

    # 경로 정의
    path = prefix
    if not os.path.isdir(path):
        os.mkdir(path)

    # mmseqs2 api 호출
    tar_gz_file = f"{path}/out.tar.gz"
    N, REDO = 101, True

    # 중복 제거 및 순서 추적
    seqs_unique = list(set(seqs))
    Ms = [N + seqs_unique.index(seq) for seq in seqs]
    # 실행!
    if not os.path.isfile(tar_gz_file):
        TIME_ESTIMATE = 150 * len(seqs_unique) # 예상 시간
        with tqdm(total=TIME_ESTIMATE, bar_format=TQDM_BAR_FORMAT) as pbar:
            while REDO:
                pbar.set_description("제출 중") # SUBMIT -> 제출 중

                # 작업이 통과될 때까지 다시 제출
                out = submit(seqs_unique, mode, N)
                while out["status"] in ["UNKNOWN", "RATELIMIT"]: # 상태가 UNKNOWN 또는 RATELIMIT인 동안
                    sleep_time = 5 + random.randint(0, 5)
                    logger.info(f"{sleep_time}초 동안 대기합니다. 이유: {out['status']}")
                    # 다시 제출
                    time.sleep(sleep_time)
                    out = submit(seqs_unique, mode, N)

                if out["status"] == "ERROR":
                    raise MMseqs2Exception()

                if out["status"] == "MAINTENANCE": # 유지보수 상태
                    raise MMseqs2Exception()

                # 작업 완료 대기
                ID, TIME = out["id"], 0
                pbar.set_description(out["status"]) # 상태 표시줄 업데이트
                while out["status"] in ["UNKNOWN", "RUNNING", "PENDING"]: # UNKNOWN, 실행 중, 대기 중 상태인 동안
                    t = 5 + random.randint(0, 5)
                    logger.info(f"{t}초 동안 대기합니다. 이유: {out['status']}")
                    time.sleep(t)
                    out = status(ID)
                    pbar.set_description(out["status"]) # 상태 표시줄 업데이트
                    if out["status"] == "RUNNING": # 실행 중인 경우
                        TIME += t
                        pbar.update(n=t)

                if out["status"] == "COMPLETE": # 완료된 경우
                    if TIME < TIME_ESTIMATE:
                        pbar.update(n=(TIME_ESTIMATE - TIME))
                    REDO = False

                if out["status"] == "ERROR": # 오류 발생 시
                    REDO = False
                    raise MMseqs2Exception()

            # 결과 다운로드
            download(ID, tar_gz_file)

    # a3m 파일 목록 준비
    if use_pairing:
        a3m_files = [f"{path}/pair.a3m"]
    else:
        a3m_files = [f"{path}/uniref.a3m"]
        if use_env:
            a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")

    # a3m 파일 추출
    if any(not os.path.isfile(a3m_file) for a3m_file in a3m_files):
        with tarfile.open(tar_gz_file) as tar_gz:
            tar_gz.extractall(path)

    # a3m 라인 수집
    a3m_lines: dict = {}
    for a3m_file in a3m_files:
        a3m_lines = get_a3m_lines(a3m_file)
    a3m_lines_list = ["".join(a3m_lines[n]) for n in Ms]

    if use_templates:
        templates_list = get_templates( # 변수명 변경: templates -> templates_list
                x,
                Path(prefix),
                "pdb70.m8",
                num_templates,
            )

    return (a3m_lines_list, templates_list) if use_templates else a3m_lines_list # 변수명 변경


def run_mmseqs_command(mmseqs: Path, params: List[Union[str, Path]]):
    module = str(params[0])
    if module in MODULE_OUTPUT_POS:
        output_pos = MODULE_OUTPUT_POS[module]
        output_path = Path(params[output_pos]).with_suffix('.dbtype')
        if output_path.exists():
            logger.info(f"{module} 모듈을 건너<0xEB><0x9B><0x84>니다. {output_path} 파일이 이미 존재합니다.")
            return

    params_log = " ".join(str(i) for i in params)
    logger.info(f"{mmseqs} {params_log} 실행 중")
    # MMseqs2의 상세 매개변수 목록이 로그를 어지럽히는 것을 숨김
    os.environ["MMSEQS_CALL_DEPTH"] = "1"
    subprocess.check_call([str(mmseqs)] + [str(i) for i in params])


# ColabFold에서 약간 수정된 코드: https://github.com/sokrypton/ColabFold
def run_local_mmseqs(
    x,
    base,
    use_env=True,
    use_templates=False,
    filter=0,
    num_templates=20,
    mmseqs_db=None,
    expand_eval: float = math.inf,
    align_eval: int = 10,
    diff: int = 3000,
    qsc: float = -20.0,
    max_accept: int = 1000000,
    prefilter_mode: int = 0,
    s: float = 8,
    db_load_mode: int = 2,
    threads: int = 32,
    gpu: int = 0,
    gpu_server: int = 0,
    unpack: bool = True,
) -> Sequence[object]:

    if filter:
        # 0.1은 위의 줄에 있는 POSIX 셸 버그로 인해 벤치마크에서 사용되지 않았습니다.
        #  EXPAND_EVAL=0.1
        align_eval = 10
        qsc = 0.8
        max_accept = 100000

    mmseqs = Path("mmseqs")
    uniref_db = Path("uniref30_2302_db")
    metagenomic_db = Path("colabfold_envdb_202108_db")
    template_db = Path("pdb100_230517")

    base.mkdir(exist_ok=True, parents=True)
    query_file = base.joinpath("query.fas")
    with query_file.open("w") as f:
        query_seq_headername = 101
        f.write(f">{query_seq_headername}\n{x}\n")

    run_mmseqs_command(
        mmseqs,
        ["createdb", query_file, base.joinpath("qdb"), "--shuffle", "0"],
    )

    used_dbs = [uniref_db]
    if use_templates:
        used_dbs.append(template_db)
    if use_env:
        used_dbs.append(metagenomic_db)

    for db in used_dbs:
        if not mmseqs_db.joinpath(f"{db}.dbtype").is_file():
            raise FileNotFoundError(f"데이터베이스 {db}가 존재하지 않습니다.")
        if (
            (
                not mmseqs_db.joinpath(f"{db}.idx").is_file()
                and not mmseqs_db.joinpath(f"{db}.idx.index").is_file()
            )
            or os.environ.get("MMSEQS_IGNORE_INDEX", False)
        ):
            logger.info("검색에서 인덱스를 사용하지 않습니다.")
            db_load_mode = 0
            dbSuffix1 = "_seq"
            dbSuffix2 = "_aln"
            dbSuffix3 = ""
        else:
            dbSuffix1 = ".idx"
            dbSuffix2 = ".idx"
            dbSuffix3 = ".idx"

    search_param = ["--num-iterations", "3",
                    "--db-load-mode", str(db_load_mode),
                    "-a", "-e", "0.1", "--max-seqs", "10000"]
    if gpu:
        # GPU 버전은 현재 갭 없는 프리필터만 지원합니다.
        search_param += ["--gpu", str(gpu), "--prefilter-mode", "1"]
    else:
        search_param += ["--prefilter-mode", str(prefilter_mode)]
        # 감도는 비 GPU 버전에 대해서만 설정할 수 있으며,
        # GPU 버전은 최대 감도로 실행됩니다.
        if s is not None:
            search_param += ["-s", "{:.1f}".format(s)]
        else:
            search_param += ["--k-score", "'seq:96,prof:80'"]
    if gpu_server:
        search_param += ["--gpu-server", str(gpu_server)]

    filter_param = ["--filter-msa", str(filter),
                    "--filter-min-enable", "1000",
                    "--diff", str(diff),
                    "--qid", "0.0,0.2,0.4,0.6,0.8,1.0",
                    "--qsc", "0", "--max-seq-id", "0.95"]
    expand_param = ["--expansion-mode", "0",
                    "-e", str(expand_eval),
                    "--expand-filter-clusters", str(filter),
                    "--max-seq-id", "0.95"]

    # uniref.a3m.dbtype이 존재하지 않으면 uniref_db 검색 실행
    if not base.joinpath("uniref.a3m").with_suffix('.a3m.dbtype').exists():
        run_mmseqs_command(mmseqs,
                           ["search", base.joinpath("qdb"),
                            mmseqs_db.joinpath(uniref_db),
                            base.joinpath("res"),
                            base.joinpath("tmp"),
                            "--threads", str(threads)] + search_param)
        run_mmseqs_command(mmseqs,
                           ["mvdb",
                            base.joinpath("tmp/latest/profile_1"),
                            base.joinpath("prof_res")])
        run_mmseqs_command(mmseqs,
                           ["lndb",
                            base.joinpath("qdb_h"),
                            base.joinpath("prof_res_h")])
        run_mmseqs_command(mmseqs,
                           ["expandaln",
                            base.joinpath("qdb"),
                            mmseqs_db.joinpath(f"{uniref_db}{dbSuffix1}"),
                            base.joinpath("res"),
                            mmseqs_db.joinpath(f"{uniref_db}{dbSuffix2}"),
                            base.joinpath("res_exp"),
                            "--db-load-mode", str(db_load_mode),
                            "--threads", str(threads)] + expand_param)
        run_mmseqs_command(mmseqs,
                           ["align",
                            base.joinpath("prof_res"),
                            mmseqs_db.joinpath(f"{uniref_db}{dbSuffix1}"),
                            base.joinpath("res_exp"),
                            base.joinpath("res_exp_realign"),
                            "--db-load-mode", str(db_load_mode),
                            "-e", str(align_eval),
                            "--max-accept", str(max_accept),
                            "--threads", str(threads),
                            "--alt-ali", "10", "-a"])
        run_mmseqs_command(mmseqs,
                           ["filterresult",
                            base.joinpath("qdb"),
                            mmseqs_db.joinpath(f"{uniref_db}{dbSuffix1}"),
                            base.joinpath("res_exp_realign"),
                            base.joinpath("res_exp_realign_filter"),
                            "--db-load-mode",
                            str(db_load_mode),
                            "--qid", "0",
                            "--qsc", str(qsc),
                            "--diff", "0",
                            "--threads", str(threads),
                            "--max-seq-id", "1.0",
                            "--filter-min-enable", "100"])
        run_mmseqs_command(mmseqs,
                           ["result2msa",
                            base.joinpath("qdb"),
                            mmseqs_db.joinpath(f"{uniref_db}{dbSuffix1}"),
                            base.joinpath("res_exp_realign_filter"),
                            base.joinpath("uniref.a3m"),
                            "--msa-format-mode", "6",
                            "--db-load-mode", str(db_load_mode),
                            "--threads", str(threads)] + filter_param)
        run_mmseqs_command(mmseqs, ["rmdb", base.joinpath("res_exp_realign_filter")])
        run_mmseqs_command(mmseqs, ["rmdb", base.joinpath("res_exp_realign")])
        run_mmseqs_command(mmseqs, ["rmdb", base.joinpath("res_exp")])
        run_mmseqs_command(mmseqs, ["rmdb", base.joinpath("res")])
    else:
        logger.info(f"uniref.a3m 파일이 이미 존재하므로 {uniref_db} 검색을 건너<0xEB><0x9B><0x84>니다.")

    # bfd.mgnify30.metaeuk30.smag30.a3m.dbtype이 존재하지 않으면 metagenomic_db 검색 실행
    bfd_exists = base.joinpath(
        "bfd.mgnify30.metaeuk30.smag30.a3m"
    ).with_suffix('.a3m.dbtype').exists()
    if use_env and not bfd_exists:
        run_mmseqs_command(mmseqs,
                           ["search",
                            base.joinpath("prof_res"),
                            mmseqs_db.joinpath(metagenomic_db),
                            base.joinpath("res_env"),
                            base.joinpath("tmp3"),
                            "--threads", str(threads)] + search_param)
        run_mmseqs_command(mmseqs,
                           ["expandaln",
                            base.joinpath("prof_res"),
                            mmseqs_db.joinpath(f"{metagenomic_db}{dbSuffix1}"),
                            base.joinpath("res_env"),
                            mmseqs_db.joinpath(f"{metagenomic_db}{dbSuffix2}"),
                            base.joinpath("res_env_exp"), "-e", str(expand_eval),
                            "--expansion-mode", "0",
                            "--db-load-mode", str(db_load_mode),
                            "--threads", str(threads)])
        run_mmseqs_command(mmseqs,
                           ["align",
                            base.joinpath("tmp3/latest/profile_1"),
                            mmseqs_db.joinpath(f"{metagenomic_db}{dbSuffix1}"),
                            base.joinpath("res_env_exp"),
                            base.joinpath("res_env_exp_realign"),
                            "--db-load-mode", str(db_load_mode),
                            "-e", str(align_eval),
                            "--max-accept", str(max_accept),
                            "--threads", str(threads),
                            "--alt-ali", "10", "-a"])
        run_mmseqs_command(mmseqs,
                           ["filterresult",
                            base.joinpath("qdb"),
                            mmseqs_db.joinpath(f"{metagenomic_db}{dbSuffix1}"),
                            base.joinpath("res_env_exp_realign"),
                            base.joinpath("res_env_exp_realign_filter"),
                            "--db-load-mode",
                            str(db_load_mode),
                            "--qid", "0",
                            "--qsc", str(qsc),
                            "--diff", "0",
                            "--max-seq-id", "1.0",
                            "--threads", str(threads),
                            "--filter-min-enable", "100"])
        run_mmseqs_command(mmseqs,
                           ["result2msa",
                            base.joinpath("qdb"),
                            mmseqs_db.joinpath(f"{metagenomic_db}{dbSuffix1}"),
                            base.joinpath("res_env_exp_realign_filter"),
                            base.joinpath("bfd.mgnify30.metaeuk30.smag30.a3m"),
                            "--msa-format-mode", "6",
                            "--db-load-mode", str(db_load_mode),
                            "--threads", str(threads)] + filter_param)
        run_mmseqs_command(mmseqs,
                           ["rmdb", base.joinpath("res_env_exp_realign_filter")])
        run_mmseqs_command(mmseqs, ["rmdb", base.joinpath("res_env_exp_realign")])
        run_mmseqs_command(mmseqs, ["rmdb", base.joinpath("res_env_exp")])
        run_mmseqs_command(mmseqs, ["rmdb", base.joinpath("res_env")])
    elif use_env:
        logger.info(
            f"bfd.mgnify30.metaeuk30.smag30.a3m 파일이 이미 존재하므로 {metagenomic_db} 검색을 건너<0xEB><0x9B><0x84>니다."
            )

    # template_db.m8.dbtype이 존재하지 않으면 template_db 검색 실행
    tmpl_db_exists = base.joinpath(
        f"{template_db}.m8"
    ).with_suffix('.m8.dbtype').exists()
    if use_templates and not tmpl_db_exists:
        run_mmseqs_command(mmseqs,
                           ["search",
                            base.joinpath("prof_res"),
                            mmseqs_db.joinpath(template_db),
                            base.joinpath("res_pdb"),
                            base.joinpath("tmp2"),
                            "--db-load-mode", str(db_load_mode),
                            "--threads", str(threads),
                            "-s", "7.5",
                            "-a", "-e", "0.1",
                            "--prefilter-mode", str(prefilter_mode)])
        run_mmseqs_command(mmseqs,
                           ["convertalis",
                            base.joinpath("prof_res"),
                            mmseqs_db.joinpath(f"{template_db}{dbSuffix3}"),
                            base.joinpath("res_pdb"),
                            base.joinpath(f"{template_db}"),
                            "--format-output",
                            "query,target,fident,alnlen,mismatch,\
gapopen,qstart,qend,tstart,tend,evalue,bits,cigar",
                            "--db-output", "1",
                            "--db-load-mode", str(db_load_mode),
                            "--threads", str(threads)])
        run_mmseqs_command(mmseqs, ["rmdb", base.joinpath("res_pdb")])
    elif use_templates:
        logger.info(
            f"{template_db}.m8 파일이 이미 존재하므로 {template_db} 검색을 건너<0xEB><0x9B><0x84>니다."
            )

    if use_env:
        run_mmseqs_command(mmseqs,
                           ["mergedbs", base.joinpath("qdb"),
                            base.joinpath("final.a3m"),
                            base.joinpath("uniref.a3m"),
                            base.joinpath("bfd.mgnify30.metaeuk30.smag30.a3m")])
        run_mmseqs_command(mmseqs,
                           ["rmdb",
                            base.joinpath("bfd.mgnify30.metaeuk30.smag30.a3m")])
        run_mmseqs_command(mmseqs,
                           ["rmdb",
                            base.joinpath("uniref.a3m")])
    else:
        run_mmseqs_command(mmseqs,
                           ["mvdb",
                            base.joinpath("uniref.a3m"),
                            base.joinpath("final.a3m")])
        run_mmseqs_command(mmseqs,
                           ["rmdb",
                            base.joinpath("uniref.a3m")])

    if unpack:
        run_mmseqs_command(mmseqs,
                           ["unpackdb",
                            base.joinpath("final.a3m"),
                            base.joinpath("."),
                            "--unpack-name-mode", "0",
                            "--unpack-suffix", ".a3m"])
        run_mmseqs_command(mmseqs,
                           ["rmdb",
                            base.joinpath("final.a3m")])

        if use_templates:
            run_mmseqs_command(mmseqs,
                               ["unpackdb",
                                base.joinpath(f"{template_db}"),
                                base.joinpath("."),
                                "--unpack-name-mode", "0",
                                "--unpack-suffix", ".m8"])
            if base.joinpath(f"{template_db}").exists():
                run_mmseqs_command(mmseqs, ["rmdb", base.joinpath(f"{template_db}")])

    run_mmseqs_command(mmseqs, ["rmdb", base.joinpath("prof_res")])
    run_mmseqs_command(mmseqs, ["rmdb", base.joinpath("prof_res_h")])
    shutil.rmtree(base.joinpath("tmp"))
    if use_templates:
        shutil.rmtree(base.joinpath("tmp2"))
    if use_env:
        shutil.rmtree(base.joinpath("tmp3"))

    if unpack:
        run_mmseqs_command(mmseqs, ["rmdb", base.joinpath("qdb")])
        run_mmseqs_command(mmseqs, ["rmdb", base.joinpath("qdb_h")])
        output_a3m = base.joinpath("0.a3m")
    else:
        output_a3m = base.joinpath("final.a3m")

    query_file.unlink()

    # API와 동일한 방식으로 a3m 라인 수집
    seqs = [x] if isinstance(x, str) else x
    N = 101
    seqs_unique = list(set(seqs))
    Ms = [N + seqs_unique.index(seq) for seq in seqs]
    a3m_lines = get_a3m_lines(output_a3m)
    a3m_lines_list = ["".join(a3m_lines[n]) for n in Ms]

    if use_templates:
        templates_list = get_templates( # 변수명 변경: templates -> templates_list
            x,
            base,
            "0.m8",
            num_templates,
            mmseqs_db=mmseqs_db,
        )

    return (a3m_lines_list, templates_list) if use_templates else a3m_lines_list # 변수명 변경


def get_a3m_lines(output_a3m):
    """
    A3M 파일을 파싱하여 서열 딕셔너리를 반환합니다.

    Args:
        output_a3m (str or Path): A3M 파일 경로.

    Returns:
        dict: 서열 ID(int)를 문자열 라인 목록(str)에 매핑하는 딕셔너리입니다.
    """
    a3m_lines: dict = {}
    update_M, M = True, None
    for line in open(output_a3m, "r"):
        if len(line) > 0:
            if "\x00" in line:
                line = line.replace("\x00", "")
                update_M = True
            if line.startswith(">") and update_M:
                M = int(line[1:].rstrip())
                update_M = False
                if M not in a3m_lines:
                    a3m_lines[M] = []
            a3m_lines[M].append(line)

    return a3m_lines


def get_templates(x, base, m8, num_templates, mmseqs_db=None):
    tested_pdbs = []
    templates_list_local = [] # 변수명 변경: templates -> templates_list_local
    logger.info("템플릿 검색 및 준비 중")
    count = 0
    for line in open(base.joinpath(m8), "r"):
        template_item = {} # 변수명 변경: template -> template_item
        if count < num_templates:
            p = line.rstrip().split()
            pdb, qid, alilen, tstart, tend = (
                p[1],
                float(p[2]),
                float(p[3]),
                int(p[8]),
                int(p[9]),
            )
            coverage = alilen / len(x)
            pdb_id = pdb.split("_")[0]

            # AF3와 동일한 템플릿 필터를 사용하고 PDB당 1개의 템플릿만 사용
            if (
                qid == 1.0
                and coverage >= 0.95
                or coverage < 0.1
                or pdb_id in tested_pdbs
            ):
                continue

            pdb_id = pdb.split("_")[0]
            if mmseqs_db:
                cif_str = fetch_local_mmcif(pdb_id,
                                            pdb.split("_")[1],
                                            tstart,
                                            tend,
                                            base,
                                            mmseqs_db)
            else:
                cif_str = fetch_mmcif(pdb_id,
                                      pdb.split("_")[1],
                                      tstart,
                                      tend,
                                      base)
            template_item["mmcif"] = cif_str # 변수명 변경

            template_seq = extract_sequence_from_mmcif(StringIO(cif_str))
            query_indices, template_indices = align_and_map(x, template_seq)

            template_item["queryIndices"] = query_indices # 변수명 변경
            template_item["templateIndices"] = template_indices # 변수명 변경
            templates_list_local.append(template_item) # 변수명 변경
            tested_pdbs.append(pdb_id)
            count += 1
    logger.info(f"다음 템플릿을 찾았습니다: {tested_pdbs}")
    return templates_list_local # 변수명 변경


def fetch_mmcif(
    pdb_id,
    chain_id,
    start,
    end,
    tmpdir,
):
    """
    주어진 PDB ID 및 체인 ID에 대한 mmCIF 파일을 가져와 AlphaFold3에서 사용할 수 있도록 준비합니다.
    """
    pdb_id = pdb_id.lower()
    url_base = "http://www.ebi.ac.uk/pdbe-srv/view/files/"
    url = url_base + pdb_id + ".cif"
    response = requests.get(url)
    text = response.text

    output = os.path.join(tmpdir, pdb_id + ".cif")
    with open(output, "w") as f:
        f.write(text)

    return get_mmcif(output, pdb_id, chain_id, start, end, tmpdir)


def fetch_local_mmcif(
        pdb_id,
        chain_id,
        start,
        end,
        tmpdir,
        mmseqs_db,
):
    """
    주어진 PDB ID 및 체인 ID에 대한 로컬 mmCIF 파일을 가져와 AlphaFold3에서 사용할 수 있도록 준비합니다.
    """
    pdb_id = pdb_id.lower()
    assert len(pdb_id) == 4, f"잘못된 PDB ID: {pdb_id}"
    inner_code = pdb_id[1:3]
    mmcif_location = mmseqs_db.joinpath(f"pdb/divided/{inner_code}/{pdb_id}.cif.gz")
    if not mmcif_location.exists():
        raise FileNotFoundError(f"MMseqs2 데이터베이스 {mmcif_location}이(가) 존재하지 않습니다.")
    with gzip.open(mmcif_location, "rb") as f:
        cif_str = f.read().decode("utf-8")

    output = os.path.join(tmpdir, pdb_id + ".cif")
    with open(output, "w") as f:
        f.write(cif_str)

    return get_mmcif(output, pdb_id, chain_id, start, end, tmpdir)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="AlphaFold3 JSON에 MMseqs2 unpaired MSA 추가"
    )
    parser.add_argument("--input_json", help="입력 AlphaFold3 JSON 파일")
    parser.add_argument("--output_json", help="출력 AlphaFold3 JSON 파일")

    parser = mmseqs2_argparse_util(parser)
    parser = custom_template_argpase_util(parser)

    args = parser.parse_args()

    add_msa_to_json(  # pragma: no cover
        args.input_json,
        args.mmseqs_database,
        args.templates,
        args.num_templates,
        False,
        args.custom_template,
        args.custom_template_chain,
        args.target_id,
        output_json=args.output_json,
        to_file=True,
    )


if __name__ == "__main__":  # pragma: no cover
    main()
