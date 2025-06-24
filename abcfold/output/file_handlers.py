import json
import logging
import re
import warnings
from abc import ABC
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from Bio.PDB import MMCIFIO, Chain, MMCIFParser, Model
from Bio.PDB.Atom import Atom
from Bio.PDB.kdtrees import KDTree
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Residue import Residue
from Bio.PDB.Superimposer import Superimposer
from Bio.SeqUtils import seq1

from abcfold.output.atoms import VANDERWALLS

warnings.filterwarnings("ignore")

logger = logging.getLogger("logger")


class FileTypes(Enum):
    """
    다양한 파일 유형에 대한 Enum 클래스입니다.
    """

    NPZ = "npz"
    NPY = "npy"
    CIF = "cif"
    JSON = "json"

    @classmethod
    def values(cls):
        return [value.value for value in cls.__members__.values()]


class ModelCount(Enum):
    """
    다양한 모델 카운트 유형에 대한 Enum 클래스입니다.
    """

    ALL = "all"
    RESIDUES = "residues"

    @classmethod
    def values(cls):
        return [value.value for value in cls.__members__.values()]


class ResidueCountType(Enum):
    """
    다양한 잔기 카운트 유형에 대한 Enum 클래스입니다.
    """

    AVERAGE = "average"
    CARBONALPHA = "carbonalpha"
    PHOSPHATE = "phosphate"

    @classmethod
    def values(cls):
        return [value.value for value in cls.__members__.values()]


class FileBase(ABC):
    """
    다양한 파일 유형에 대한 추상 기본 클래스입니다.
    """

    def __init__(self, pathway: Union[str, Path]):
        self.pathway = Path(pathway)
        self.suffix = self.pathway.suffix[1:]

    def __str__(self):
        return str(self.pathway)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.pathway})"


class NpzFile(FileBase):

    def __init__(self, npz_file: Union[str, Path]):
        """
        npz 파일을 처리하는 객체입니다.

        Args:
            npz_file (Union[str, Path]): npz 파일 경로

        Attributes:
            npz_file (Path): npz 파일 경로
            data (dict): npz 파일의 데이터를 포함하는 딕셔너리

        """
        super().__init__(npz_file)
        self.npz_file = Path(npz_file)
        self.data = self.load_npz_file()

    def load_npz_file(self) -> dict:
        return dict(np.load(self.npz_file))


class NpyFile(FileBase):
    def __init__(self, npy_file: Union[str, Path]):
        """
        npy 파일을 처리하는 객체입니다.

        Args:
            npy_file (Union[str, Path]): npy 파일 경로

        Attributes:
            npy_file (Path): npy 파일 경로
            data (np.ndarray): npy 파일의 데이터를 포함하는 Numpy 배열
        """

        super().__init__(npy_file)
        self.npy_file = Path(npy_file)
        self.data = self.load_npy_file()

    def load_npy_file(self) -> np.ndarray:
        return np.load(self.npy_file)


class CifFile(FileBase):

    def __init__(self, cif_file: Union[str, Path], input_params: Optional[dict] = None):
        """
        cif 파일을 처리하는 객체입니다.

        Args:
            cif_file (Union[str, Path]): cif 파일 경로
            input_params (Optional[dict]): 모델에 사용된 입력 매개변수를 포함하는 딕셔너리입니다.
                                         리간드와 서열을 구별하는 데 사용됩니다.

        Attributes:
            cif_file (Path): cif 파일 경로
            model (Structure): 모델을 포함하는 BioPython 구조 객체
            atom_plddt_per_chain (dict): 각 원자에 대한 pLDDT 점수를 포함하는 딕셔너리
            residue_plddt_per_chain (dict): 각 잔기에 대한 pLDDT 점수를 포함하는 딕셔너리
            plddts (list): 각 원자에 대한 pLDDT 점수를 포함하는 목록
            residue_plddts (list): 각 잔기에 대한 pLDDT 점수를 포함하는 목록
            name (str): 모델에 지정된 이름
        """
        if input_params is None:
            self.input_params = {}
        else:
            self.input_params = input_params

        super().__init__(cif_file)
        self.cif_file = Path(cif_file)
        self.clashes = 0
        self.clashes_residues = 0
        self.model = self.load_cif_file()
        self.__ligand_plddts = None
        self.__plddts = None
        self.__residue_plddts = None
        self.__h_score = None
        self.__name = self.cif_file.stem

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, name: str):
        if not isinstance(name, str):
            logger.error("이름은 문자열이어야 합니다.")
            raise ValueError()
        self.__name = name

    @property
    def plddts(self):
        """
        모델의 각 원자에 대한 pLDDT 점수입니다.
        """
        self.__plddts = [
            plddts for plddts in self.get_plddt_per_atom().values() for plddts in plddts
        ]
        return self.__plddts

    @property
    def residue_plddts(self):
        """
        모델의 각 잔기에 대한 pLDDT 점수입니다.
        """

        self.__residue_plddts = [
            plddts
            for plddts in self.get_plddt_per_residue().values()
            for plddts in plddts
        ]
        return self.__residue_plddts

    @property
    def average_plddt(self):
        """
        모델의 평균 pLDDT 점수입니다.
        """
        return float(np.mean(self.plddts))

    @property
    def ligand_plddts(self):
        """
        모델의 각 리간드에 대한 pLDDT 점수입니다.
        """
        self.__ligand_plddts = self.get_plddt_per_ligand()
        return self.__ligand_plddts

    @property
    def h_score(self):
        """
        모델의 H 점수입니다.
        """
        self.__h_score = self.calculate_h_score()
        return self.__h_score

    def load_cif_file(self):
        """
        BioPython을 사용하여 cif 파일을 로드합니다.
        """
        parser = MMCIFParser(QUIET=True)
        return parser.get_structure(self.pathway.stem, self.pathway)

    def get_chains(self):
        """
        모델의 체인 목록을 가져옵니다.
        """
        return self.model[0]

    def chain_lengths(
        self,
        mode=ModelCount.RESIDUES,
        ligand_atoms=False,
        ptm_atoms=False,
    ) -> dict:
        """
        모델의 각 체인 길이를 가져오는 함수입니다.

        Args:
            mode (ModelCount): 사용할 모드를 지정하는 Enum 클래스입니다.
                               참고: 리간드의 경우 길이는 항상 원자 수가 됩니다.

        Returns:
            dict: 체인 ID와 체인 길이를 포함하는 딕셔너리입니다.

        Raises:
            ValueError: 모드가 유효하지 않은 경우 발생합니다.
        """
        chains = self.get_chains()
        if mode == ModelCount.ALL or mode == ModelCount.ALL.value:

            # return {
            #     chain.id: len([atom for resiude in chain for atom in resiude])
            #     for chain in chains
            # }
            chain_lengths: dict = {}
            for chain in chains:
                if chain.id in chain_lengths:
                    chain_lengths[chain.id] += len(
                        [atom for resiude in chain for atom in resiude]
                    )
                else:
                    chain_lengths[chain.id] = len(
                        [atom for resiude in chain for atom in resiude]
                    )

            return chain_lengths

        elif mode == ModelCount.RESIDUES or mode == ModelCount.RESIDUES.value:
            residue_counts: dict = {}
            for chain in chains:
                if self.check_other(chain, ["protein", "rna", "dna"]) and ptm_atoms:
                    counter = 0
                    for residue in chain:
                        if is_aa(residue.resname, standard=True):
                            counter += 1
                            continue
                        else:
                            counter += len([atom for atom in residue])

                    residue_counts[chain.id] = counter

                elif self.check_ligand(chain):
                    if ligand_atoms:
                        if chain.id in residue_counts:
                            residue_counts[chain.id] += len(
                                [atom for resiude in chain for atom in resiude]
                            )
                        else:
                            residue_counts[chain.id] = len(
                                [atom for resiude in chain for atom in resiude]
                            )

                    else:
                        residue_counts[chain.id] = 1
                else:
                    residue_counts[chain.id] = len([residue for residue in chain])

            return residue_counts

        else:
            msg = f"유효하지 않은 모드입니다. {', '.join(ModelCount.__members__)} 중 하나를 사용하십시오."
            logger.critical(msg)
            raise ValueError()

    def token_residue_ids(self) -> dict:
        """
        모델의 각 체인에 대한 잔기 ID를 가져오는 함수입니다.

        Returns:
            dict: 각 체인의 체인 ID와 잔기 ID를 포함하는 딕셔너리입니다.
        """
        from abcfold.output.utils import flatten

        chains = self.get_chains()
        residue_ids = {}
        for chain in chains:
            if self.check_ligand(chain):
                residue_ids[chain.id] = [
                    [residue.id[1]] for residue in chain for _ in residue
                ]
                continue
            residue_ids[chain.id] = [
                (
                    [residue.id[1]]
                    if residue.id[0] == " " or residue.id[0] == "H"
                    else [residue.id[1] for _ in residue]
                )
                for residue in chain
            ]

        residue_ids = {k: flatten(v) for k, v in residue_ids.items()}

        return residue_ids

    def calculate_h_score(self):
        """
        모델의 H 점수를 계산합니다.

        Returns:
            float: 모델의 H 점수입니다.
        """

        score = 0
        for i in reversed(range(1, 101)):
            if (100.0 / len(self.plddts)) * np.sum(np.array(self.plddts) >= i) >= i:
                score = i
                break
        return score

    def get_model_sequence_data(self) -> dict:
        """
        모델의 각 체인 및 리간드에 대한 서열을 가져옵니다. 플로팅을 위해 내부적으로 사용됩니다.

        Returns:
            dict : 체인 ID 및 서열 데이터
        """
        sequence_data = {}
        for chain in self.model[0]:
            if self.check_ligand(chain):
                sequence_data[chain.id] = "".join(
                    [atom.id[0] for residue in chain for atom in residue]
                )
            else:
                sequence_data[chain.id] = "".join(
                    [seq1(residue.get_resname()) for residue in chain]
                )
        return sequence_data

    def get_plddt_per_atom(self) -> dict:
        """
        모델의 각 원자에 대한 pLDDT 점수를 가져옵니다.

        Returns:
            dict: 각 원자의 체인 ID와 pLDDT 점수를 포함하는 딕셔너리입니다.
        """
        plddt: Dict[str, list] = {}
        for chain in self.model[0]:

            for residue in chain:
                for atom in residue:
                    if chain.id in plddt:
                        plddt[chain.id].append(atom.bfactor)
                    else:
                        plddt[chain.id] = [atom.bfactor]

        return plddt

    def get_plddt_per_residue(self, method=ResidueCountType.AVERAGE.value) -> dict:
        """
        모델의 각 잔기에 대한 pLDDT 점수를 가져옵니다.

        Args:
            method (ResidueCountType): 사용할 방법을 지정하는 Enum 클래스입니다.

        Returns:
            dict: 각 잔기의 체인 ID와 pLDDT 점수를 포함하는 딕셔너리입니다.
        """
        plddts: Dict[str, list] = {}

        if method not in ResidueCountType.values():
            logger.error(
                f"유효하지 않은 방법입니다. {', '.join(ResidueCountType.__members__)} 중 하나를 사용하십시오."
            )
            raise ValueError()

        chains = self.get_chains()
        for chain in chains:
            if self.check_ligand(chain):
                if chain.id in plddts:
                    plddts[chain.id].extend(
                        [atom.bfactor for residue in chain for atom in residue]
                    )
                else:
                    plddts[chain.id] = [
                        atom.bfactor for residue in chain for atom in residue
                    ]

            else:
                for residue in chain:
                    if method == ResidueCountType.AVERAGE.value:
                        scores = 0
                        for atom in residue:
                            scores += atom.bfactor
                        score = scores / len(residue)

                    elif method == ResidueCountType.CARBONALPHA.value:
                        for atom in residue:
                            if atom.id == "CA":
                                score = atom.bfactor
                                break

                    elif method == ResidueCountType.PHOSPHATE.value:
                        for atom in residue:
                            if atom.id == "P":
                                score = atom.bfactor
                                break

                    if chain.id in plddts:
                        plddts[chain.id].append(score)

                    else:
                        plddts[chain.id] = [score]

        plddt_lengths = {k: len(v) for (k, v) in plddts.items()}
        chain_lengths = self.chain_lengths(mode="residues", ligand_atoms=True)
        for chain_id in plddt_lengths:
            assert (
                chain_lengths[chain_id] == plddt_lengths[chain_id]
            ), f"{chain_id}, {chain_lengths[chain_id]} != {plddt_lengths[chain_id]}" # 길이 불일치 확인
        return plddts

    def get_plddt_per_ligand(self) -> dict:
        """
        모델의 각 리간드에 대한 pLDDT 점수를 가져옵니다.

        Returns:
            dict: 각 원자의 체인 ID와 pLDDT 점수를 포함하는 딕셔너리입니다.
        """
        plddt: Dict[str, list] = {}
        for chain in self.model[0]:
            if self.check_ligand(chain):
                for residue in chain:
                    for atom in residue:
                        if chain.id in plddt:
                            plddt[chain.id].append(atom.bfactor)
                        else:
                            plddt[chain.id] = [atom.bfactor]
        return plddt

    def check_ligand(self, chain: Chain) -> bool:
        """
        체인이 리간드인지 확인합니다.

        Args:
            chain (Chain): BioPython 체인 객체

        Returns:
            bool: 체인이 리간드이면 True, 그렇지 않으면 False
        """

        return self.check_other(chain, ["ligand"])

    def check_other(self, chain: Chain, check_list) -> bool:
        """
        체인이 check_list의 유형 중 하나인지 확인합니다.

        Args:
            chain (Chain): 확인할 BioPython 체인 객체
            check_list (list): 확인할 유형 목록 (예: ["protein", "dna"])

        Returns:
            bool: 체인이 check_list의 유형 중 하나이면 True, 그렇지 않으면 False
        """
        sequences = self.input_params.get("sequences")
        if sequences is None:
            logger.warning("입력 파일에서 서열 정보를 가져올 수 없습니다.")
            return False
        for sequence in sequences:
            for sequence_type, sequence_data in sequence.items():
                if sequence_type in check_list:
                    if "id" not in sequence_data:
                        continue
                    if hasattr(chain, "id"):
                        chain_id = chain.id
                    else:
                        chain_id = chain
                    if isinstance(sequence_data["id"], str):
                        if chain_id == sequence_data["id"]:
                            return True
                    elif isinstance(sequence_data["id"], list):
                        if chain_id in sequence_data["id"]:
                            return True
        return False

    def relabel_chains(
        self, chain_ids: List[str], link_ids: Optional[dict] = None
    ) -> None:
        """
        모델의 체인 레이블을 변경합니다.

        Args:
            chain_ids (List[str]): 체인 레이블을 변경할 체인 ID 목록입니다.
            link_ids (Optional[dict]): 연결된 ID들을 나타내는 딕셔너리입니다.

        Returns:
            None
        """

        chain_ids = chain_ids.copy()
        structure = self.model[0]
        old_new_chain_id = {}

        if link_ids is None:
            link_ids = {}
        else:
            for new_ids in link_ids.values():
                for new_id in new_ids:
                    chain_ids.pop(chain_ids.index(new_id))

        old_chain_label_counter, new_chain_label_counter = 0, 0

        chain_names = [chain.id for chain in structure]
        while old_chain_label_counter < len(structure):
            chain = chain_ids[new_chain_label_counter]
            old_new_chain_id[chain_names[old_chain_label_counter]] = chain

            # 체인이 레이블 변경될 때마다 old_chain 증가
            old_chain_label_counter += 1

            if chain in link_ids:
                ligand_no_added = 2
                for _ in link_ids[chain]:
                    old_new_chain_id[chain_names[old_chain_label_counter]] = chain
                    for residue in structure[chain_names[old_chain_label_counter]]:

                        residue.id = (residue.id[0], ligand_no_added, residue.id[2])
                        ligand_no_added += 1
                    old_chain_label_counter += 1

            new_chain_label_counter += 1

        for chain_to_rename in structure:
            chain_to_rename.id = old_new_chain_id[chain_to_rename.id]

        assert old_chain_label_counter == len(
            self.get_chains()
        ), "체인 ID 수는 체인 수와 일치해야 합니다."
        self.update()

    def update(self):
        """
        현재 CifFile 객체를 파일에 쓰고 다시 로드하여 업데이트합니다.
        """
        self.to_file(self.pathway)
        self = CifFile(self.pathway, self.input_params)

    def reorder_chains(self, new_chain_ids: List[str]):
        """
        모델의 체인 순서를 변경합니다.

        Args:
            new_chain_ids (List[str]): 새로운 체인 ID 순서 목록입니다.
        """

        assert sorted([chain.id for chain in self.get_chains()]) == sorted(
            new_chain_ids
        ), "재정렬을 위해서는 체인 ID가 모델에 이미 있는 것과 동일해야 합니다."

        new_model = Model.Model(0)

        [new_model.add(self.model[0][ch]) for ch in new_chain_ids]
        self.model.detach_child(0)
        self.model.add(new_model)
        self.update()

    def check_clashes(
        self,
        threshold: Union[int, float] = 3.4,
        bucket: int = 10,
        clash_cutoff: float = 0.63,
    ) -> Tuple[List[Tuple[Atom, Atom]], List[Tuple[Residue, Residue]]]:
        """
        다른 체인의 원자 간 충돌을 확인합니다.

        Args:
            threshold: 충돌에 대한 거리 임계값입니다.
            bucket: KDTree 버킷 크기입니다.
            clash_cutoff: 반 데르 발스 반경에 대한 충돌 계수입니다.

        Returns:
            충돌 목록입니다. (원자 쌍 목록, 잔기 쌍 목록)

        """
        atoms = self.get_atoms()
        coords = np.array(
            [atom.get_coord() for atom in atoms],
            dtype="d",
        )
        assert bucket > 1
        assert coords.shape[1] == 3
        assert clash_cutoff > 0.0 and clash_cutoff <= 1.0

        tree = KDTree(coords, bucket)
        neighbors = tree.neighbor_search(threshold)
        clashes_atoms, clashes_residues = [], []

        for neighbor in neighbors:
            i1, i2 = neighbor.index1, neighbor.index2
            atom1, atom2 = atoms[i1], atoms[i2]
            # 원자의 원소 가져오기
            element1 = atom1.element
            element2 = atom2.element
            # chain_id 및 residue_id 찾기
            chain_id1 = atom1.get_full_id()[2]
            chain_id2 = atom2.get_full_id()[2]

            if chain_id1 == chain_id2:
                continue

            distance = np.linalg.norm(atom1.get_coord() - atom2.get_coord())
            if (atom1.name == "C" and atom1.name == "N") or (
                atom2.name == "N" and atom1.name == "C"
            ):
                continue
            elif (atom1.name == "SG" and atom2.name == "SG") and distance > 1.88:
                continue

            clash_radius = (
                VANDERWALLS.get(element1, 1.7) + VANDERWALLS.get(element2, 1.7)
            ) * 0.63
            if distance < clash_radius:
                residue1 = atom1.get_parent()
                residue2 = atom2.get_parent()

                clashes_atoms.append((atom1, atom2))

                if (residue1, residue2) not in clashes_residues:
                    clashes_residues.append((residue1, residue2))

        self.clashes = len(clashes_atoms)
        self.clashes_residues = len(clashes_residues)
        return (clashes_atoms, clashes_residues)

    def get_atoms(self, chain_id=None) -> list:
        """
        구조의 원자를 가져옵니다.

        Args:
            chain_id (str, optional): 특정 체인의 원자만 가져오려면 체인 ID를 지정합니다. 기본값은 None입니다.

        Returns:
            list: 원자 목록입니다.

        """
        if chain_id is not None:
            return [
                atom
                for chain in self.model[0]
                for atom in chain.get_atoms()
                if chain.id == chain_id
            ]
        return [atom for chain in self.model[0] for atom in chain.get_atoms()]

    def to_file(self, output_file: Union[str, Path]) -> None:
        """
        cif 파일을 저장합니다.

        Args:
            output_file (Union[str, Path]): cif 파일을 저장할 경로입니다.

        Returns:
            None
        """
        io = MMCIFIO()
        io.set_structure(self.model)

        # save는 딕셔너리를 생성합니다.
        io.save(str(output_file))
        self.__atom_site_label_update(io.dic)
        self.__ligand_to_hetatm(io.dic)
        with open(output_file, "w") as f:
            io._save_dict(f)

        self.__single_to_double_quotes(output_file)

    def __single_to_double_quotes(self, file_name: Union[str, Path]) -> None:
        # 파일 내 작은따옴표로 둘러싸인 문자열을 큰따옴표로 변경 (특정 PDB 형식 문제 해결용)
        new_lines = []
        with open(file_name, "r") as f:
            lines = [line.rstrip() for line in f]

        for line in lines:

            single_quotes = re.compile(r"[']\w+[']{2}")
            new_lines.append(
                single_quotes.sub(lambda x: f'"{x.group()[1:-1]}"', line, count=1)
            )

        with open(file_name, "w") as f:
            f.write("\n".join(new_lines))

    def __atom_site_label_update(self, out_dict):
        # _atom_site.label_asym_id 필드를 체인 길이에 맞게 업데이트
        atom_site_labels_asym_ids = []
        for chain_id, chain_length in self.chain_lengths(mode="all").items():
            atom_site_labels_asym_ids.extend([chain_id] * chain_length)

        assert len(out_dict["_atom_site.label_asym_id"]) == len(
            atom_site_labels_asym_ids
        ), f"Lengths must be the same, current lengths are \
{len(out_dict['_atom_site.label_asym_id'])} and {len(atom_site_labels_asym_ids)}"

        out_dict["_atom_site.label_asym_id"] = atom_site_labels_asym_ids

        return out_dict

    def __ligand_to_hetatm(self, out_dict):
        # 리간드 체인의 _atom_site.group_PDB 필드를 "HETATM"으로 설정
        atom_site_group_pdb = []
        counter = 0
        for chain_id, chain_length in self.chain_lengths(mode="all").items():
            if self.check_ligand(chain_id):
                atom_site_group_pdb.extend(["HETATM"] * chain_length)
            else:
                atom_site_group_pdb.extend(
                    out_dict["_atom_site.group_PDB"][
                        counter : counter + chain_length  # noqa: E203
                    ]
                )

            counter += chain_length

        assert len(out_dict["_atom_site.group_PDB"]) == len(atom_site_group_pdb)

        out_dict["_atom_site.group_PDB"] = atom_site_group_pdb

        return out_dict


class ConfidenceJsonFile(FileBase):
    def __init__(self, json_file: Union[str, Path]):
        """
        JSON 파일을 처리하는 객체입니다. (주로 신뢰도 점수 관련)

        Args:
            json_file (Union[str, Path]): JSON 파일 경로

        Attributes:
            json_file (Path): JSON 파일 경로
            data (dict): JSON 파일의 데이터를 포함하는 딕셔너리

        """
        super().__init__(json_file)
        self.data = self.load_json_file()

    def load_json_file(self):
        # JSON 파일 로드
        with open(self.pathway, "r") as f:
            data = json.load(f)

        return data


def superpose_models(models_list: List[Union[str, Path]]) -> None:
    """
    목록의 모델들을 중첩시키고 새 파일에 저장합니다.

    Args:
        models_list (List[Union[str, Path]]): 중첩할 모델 목록

    Returns:
        None
    """

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(Path(models_list[0]).stem, Path(models_list[0]))
    ref_model = structure[0]

    for model in models_list[1:]:
        alt_structure = parser.get_structure(Path(model).stem, Path(model))
        alt_model = alt_structure[0]

        ref_atoms = []
        alt_atoms = []
        for (ref_chain, alt_chain) in zip(ref_model, alt_model):
            for ref_res, alt_res in zip(ref_chain, alt_chain):
                if ref_res.resname != alt_res.resname or ref_res.id != alt_res.id:
                    pass

                # 뉴클레오타이드와 단백질을 다르게 처리
                if ref_res.resname in ["DA", "DT", "DG", "DC"]:
                    ref_atoms.append(ref_res["C1'"])
                    alt_atoms.append(alt_res["C1'"])
                elif ref_res.resname in ["A", "U", "G", "C", "T"]:
                    ref_atoms.append(ref_res["C1'"])
                    alt_atoms.append(alt_res["C1'"])
                elif 'CA' in ref_res:
                    ref_atoms.append(ref_res['CA'])
                    alt_atoms.append(alt_res['CA'])
                else:  # 그 외 다른 것은 무시
                    pass

        if len(ref_atoms) == 0 or len(alt_atoms) == 0:
            logger.warning(
                f"{model}에서 중첩을 위한 일치하는 원자를 찾을 수 없습니다. 건너<0xEB><0x9B><0x84>니다."
            )
        else:
            super_imposer = Superimposer()
            super_imposer.set_atoms(ref_atoms, alt_atoms)
            super_imposer.apply(alt_model.get_atoms())

            io = MMCIFIO()
            io.set_structure(alt_structure)
            io.save(str(model))  # 원본 파일 덮어쓰기
