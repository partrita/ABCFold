import json
from itertools import zip_longest
from pathlib import Path
from typing import List, Union

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner

from abcfold.output.file_handlers import CifFile

AF3TEMPLATE: dict = {
    "atom_chain_ids": [],
    "atom_plddts": [],
    "contact_probs": [],
    "pae": [],
    "token_chain_ids": [],
    "token_res_ids": [],
}


class Af3Pae:
    @classmethod
    def from_alphafold3(cls, scores: dict, cif_file: CifFile):
        """AlphaFold3 출력에서 Af3Pae 객체를 생성합니다."""
        def reorder_matrix(pae_matrix, chain_lengths, af3_chain_lengths):
            # PAE 행렬이 numpy 배열이 아니면 변환합니다.
            if not isinstance(pae_matrix, np.ndarray):
                pae_matrix = np.array(pae_matrix)
            desired_order = flatten(
                [[key] * value for key, value in chain_lengths.items()]
            )
            current_order = flatten(
                [[key] * value for key, value in af3_chain_lengths.items()]
            )

            order_mapping = {}

            for i, chain_id in enumerate(current_order):
                if chain_id in desired_order:
                    order_mapping[i] = desired_order.index(chain_id)
                    desired_order[desired_order.index(chain_id)] = None  # Mark as used
                else:
                    order_mapping[i] = i

            # Reorder the PAE matrix rows and columns based on the mapping

            reordered_matrix = np.zeros_like(pae_matrix)
            for i in range(len(pae_matrix)):
                for j in range(len(pae_matrix)):

                    reordered_matrix[order_mapping[i], order_mapping[j]] = pae_matrix[
                        i, j
                    ]

            return reordered_matrix.tolist()

        af3_scores = AF3TEMPLATE.copy()

        chain_lengths = cif_file.chain_lengths(
            mode="residues", ligand_atoms=True, ptm_atoms=True
        )

        af3pae_chain_lengths = {
            k: chain_lengths[k] for k in np.unique(scores["token_chain_ids"])
        }
        if list(chain_lengths.keys()) == list(af3pae_chain_lengths.keys()):
            return cls(scores)

        residue_lengths = cif_file.chain_lengths(mode="all", ligand_atoms=True)

        atom_chain_ids = flatten(
            [[key] * value for key, value in residue_lengths.items()]
        )
        atom_plddts = cif_file.plddts

        token_res_ids_dict = cif_file.token_residue_ids()
        token_res_ids = flatten([value for _, value in token_res_ids_dict.items()])

        reordered_pae = reorder_matrix(
            scores["pae"], chain_lengths, af3pae_chain_lengths
        )
        contact_probs = reorder_matrix(
            scores["contact_probs"], chain_lengths, af3pae_chain_lengths
        )
        token_chain_ids = flatten(
            [[key] * len(value) for key, value in token_res_ids_dict.items()]
        )

        assert len(atom_chain_ids) == len(scores["atom_chain_ids"])
        assert sum([len(value) for value in token_res_ids_dict.values()]) == len(
            scores["pae"]
        )
        assert len(atom_plddts) == len(scores["atom_plddts"])

        assert len(reordered_pae) == len(scores["pae"])
        assert len(contact_probs) == len(scores["contact_probs"])
        assert len(reordered_pae[0]) == len(scores["pae"][0])
        assert len(contact_probs[0]) == len(scores["contact_probs"][0])
        assert len(token_chain_ids) == len(scores["token_chain_ids"])
        assert len(token_res_ids) == len(scores["token_res_ids"])

        af3_scores["atom_chain_ids"] = atom_chain_ids
        af3_scores["atom_plddts"] = atom_plddts
        af3_scores["contact_probs"] = contact_probs
        af3_scores["pae"] = reordered_pae
        af3_scores["token_chain_ids"] = token_chain_ids
        af3_scores["token_res_ids"] = token_res_ids

        return cls(af3_scores)

    @classmethod
    def from_boltz(cls, scores: dict, cif_file: CifFile):
        """Boltz 출력에서 Af3Pae 객체를 생성합니다."""
        af3_scores = AF3TEMPLATE.copy()

        chain_lengths = cif_file.chain_lengths(mode="residues", ligand_atoms=True)
        residue_lengths = cif_file.chain_lengths(mode="all", ligand_atoms=True)

        atom_chain_ids = flatten(
            [[key] * value for key, value in residue_lengths.items()]
        )

        atom_plddts = cif_file.plddts
        token_chain_ids = flatten(
            [[key] * value for key, value in chain_lengths.items()]
        )

        token_res_ids = flatten(
            [
                [value for value in values]
                for _, values in cif_file.token_residue_ids().items()
            ]
        )

        af3_scores["pae"] = scores["pae"].tolist()
        af3_scores["atom_chain_ids"] = atom_chain_ids
        af3_scores["atom_plddts"] = atom_plddts
        af3_scores["contact_probs"] = np.zeros(shape=scores["pae"].shape).tolist()
        af3_scores["token_chain_ids"] = token_chain_ids
        af3_scores["token_res_ids"] = token_res_ids

        return cls(af3_scores)

    @classmethod
    def from_chai1(cls, scores: np.ndarray, cif_file: CifFile):
        """Chai-1 출력에서 Af3Pae 객체를 생성합니다."""
        af3_scores = AF3TEMPLATE.copy()
        chain_lengths = cif_file.chain_lengths(mode="residues", ligand_atoms=True)

        residue_lengths = cif_file.chain_lengths(mode="all", ligand_atoms=True)

        atom_chain_ids = flatten(
            [[key] * value for key, value in residue_lengths.items()]
        )

        atom_plddts = cif_file.plddts
        token_chain_ids = flatten(
            [[key] * value for key, value in chain_lengths.items()]
        )

        token_res_ids = flatten(
            [
                [value for value in values]
                for _, values in cif_file.token_residue_ids().items()
            ]
        )

        af3_scores["pae"] = scores.tolist()
        af3_scores["atom_chain_ids"] = atom_chain_ids
        af3_scores["atom_plddts"] = atom_plddts
        af3_scores["contact_probs"] = np.zeros(shape=scores.shape).tolist()
        af3_scores["token_chain_ids"] = token_chain_ids
        af3_scores["token_res_ids"] = token_res_ids

        return cls(af3_scores)

    def __init__(self, af3_scores: dict):
        """Af3Pae 객체를 초기화합니다."""
        self.scores = af3_scores

    def to_file(self, file_path: Union[str, Path]):
        """PAE 점수를 JSON 파일로 저장합니다."""
        with open(file_path, "w") as f:
            json.dump(self.scores, f, indent=4)


def flatten(xss):
    """중첩된 리스트를 단일 리스트로 평탄화합니다."""
    return [x for xs in xss for x in xs]


def get_gap_indicies(*cif_objs) -> List[np.ndarray]:
    """
    CIF 객체들 사이의 갭을 가져옵니다. 모델링 프로그램 간에 체인 길이가
    일치하지 않는 경우가 있는데, 이 함수는 이러한 불일치가 발생하는 위치를 찾는 데 사용됩니다.

    Args:
        *cif_objs: 다중 CIF 객체

    Returns:
        indicies: 불일치가 위치한 체인 ID를 키로 하고, -1이 갭을 나타내는
                  인덱스 목록을 값으로 하는 딕셔너리입니다.
                  (현재 반환 타입은 List[np.ndarray]로 되어 있으나, 설명은 딕셔너리 형태로 되어 있어 확인 필요)
    """
    indicies: list = []

    if len(cif_objs) == 1:
        return indicies
    chain_lengths = [
        cif.chain_lengths(mode="residues", ligand_atoms=True) for cif in cif_objs
    ]

    assert all(
        [
            chain_lengths[0].keys() == chain_lengths[i].keys()
            for i in range(1, len(chain_lengths) - 1)
        ]
    )

    unequal_chain_lengths = [
        id_
        for id_ in chain_lengths[0].keys()
        if any(
            [
                chain_lengths[0][id_] != chain_lengths[i][id_]
                for (i, _) in enumerate(chain_lengths[1:], start=1)
            ]
        )
    ]

    for chain_id in chain_lengths[0]:
        if chain_id in unequal_chain_lengths:
            chain_atoms = [
                "".join([atom.element for atom in cif.get_atoms(chain_id=chain_id)])
                for cif in cif_objs
            ]

            longest = max(chain_atoms, key=len)

            for atom_str in chain_atoms:
                alignment = PairwiseAligner().align(longest, atom_str)
                indicies.append(alignment[0].indices[1])
        else:
            for _ in cif_objs:

                indicies.append(np.array([1] * chain_lengths[0][chain_id]))

    indicies = interleave_repeated(
        indicies, len(cif_objs), len(list(chain_lengths[0].keys()))
    )

    return indicies


def interleave_repeated(lst, n, chain_no):
    """
    반복을 기반으로 리스트의 리스트(청크)에서 요소들을 인터리빙합니다.
    get_gap_indicies의 헬퍼 함수입니다.

    Args:
        lst (list): 처리될 리스트로, 청크에서 평탄화된 것으로 예상됩니다.
        n (int): 원본 CIF 객체의 수 (반복 횟수).
        chain_no (int): 체인의 수.

    Returns:
        list: 원본 CIF 객체별로 그룹화된 인터리빙된 인덱스 목록입니다.
    """
    indicies = []
    chunks = [lst[i : i + n] for i in range(0, len(lst), n)]  # noqa: E203
    interleaved = [x for tup in zip_longest(*chunks) for x in tup if x is not None]

    for i in range(0, len(interleaved), chain_no):
        tmp_lst = []
        for j in range(chain_no):
            tmp_lst.extend(interleaved[i + j])
        indicies.append(tmp_lst)

    return indicies


def insert_none_by_minus_one(indices, values):
    """
    `indices`의 해당 인덱스가 -1인 모든 위치에 None을 값 목록에 삽입합니다.

    Args:
        indices (list): -1이 갭을 나타내는 인덱스 목록입니다.
        values (list): 삽입될 값 목록입니다.

    Returns:
        list: 갭 위치에 None이 삽입된 새 목록입니다.
    """
    result = []
    value_index = 0

    for idx in indices:
        if idx == -1:
            result.append(None)
        else:
            result.append(values[value_index])
            value_index += 1

    assert len(indices) == len(result)

    return result


def make_dummy_m8_file(run_json, output_dir):
    """
    실행 JSON 파일의 템플릿으로 더미 m8 파일을 만듭니다.
    """
    with open(run_json) as f:
        input_json = json.load(f)

    templates = {}
    for sequence in input_json["sequences"]:
        if "protein" not in sequence:
            continue
        for id_ in sequence["protein"]["id"]:
            for template in sequence["protein"]["templates"]:
                if id_ in templates:
                    templates[id_].append(
                        template["mmcif"].split("\n")[0].split("_")[1]
                    )
                else:
                    templates[id_] = [template["mmcif"].split("\n")[0].split("_")[1]]

    m8_file = output_dir / "dummy.m8"
    if not templates:
        return None
    table = []

    for id_ in templates:
        for template in templates[id_]:
            table.append([id_, template, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    pd.DataFrame(table).to_csv(m8_file, sep="\t", header=False, index=False)

    return m8_file
