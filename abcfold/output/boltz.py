import logging
from pathlib import Path
from typing import Union

from abcfold.boltz.af3_to_boltz import BoltzYaml
from abcfold.output.file_handlers import (CifFile, ConfidenceJsonFile,
                                          FileTypes, ModelCount, NpzFile)
from abcfold.output.utils import Af3Pae

logger = logging.getLogger("logger")


class BoltzOutput:
    def __init__(
        self,
        boltz_output_dir: Union[str, Path],
        input_params: dict,
        name: str,
    ):
        """
        Boltz 실행 출력을 처리하는 객체입니다.

        Args:
            boltz_output_dir (Union[str, Path]): Boltz 출력 디렉터리 경로
            input_params (dict): Boltz 실행에 사용된 입력 매개변수를 포함하는 딕셔너리
            name (str): Boltz 실행에 지정된 이름

        Attributes:
            output_dir (Path): Boltz 출력 디렉터리 경로
            input_params (dict): Boltz 실행에 사용된 입력 매개변수를 포함하는 딕셔너리
            name (str): Boltz 실행에 지정된 이름
            output (dict): Boltz 출력 디렉터리의 내용을 처리한 출력을 포함하는 딕셔너리입니다.
                           딕셔너리 구조는 다음과 같습니다:

            {
                1: {
                    "pae": NpzFile,
                    "plddt": NpzFile,
                    "pde": NpzFile,
                    "cif": CifFile,
                    "json": ConfidenceJsonFile
                },
                2: {
                    "pae": NpzFile,
                    "plddt": NpzFile,
                    "pde": NpzFile,
                    "cif": CifFile,
                    "json": ConfidenceJsonFile
                },
                ...
            }
            pae_files (list): PAE 데이터를 포함하는 정렬된 NpzFile 객체 목록
            plddt_files (list): PLDDT 데이터를 포함하는 정렬된 NpzFile 객체 목록
            pde_files (list): PDE 데이터를 포함하는 정렬된 NpzFile 객체 목록
            cif_files (list): 모델 데이터를 포함하는 정렬된 CifFile 객체 목록
            scores_files (list): 모델 점수를 포함하는 정렬된 ConfidenceJsonFile 객체 목록

        """
        self.output_dir = Path(boltz_output_dir)
        self.input_params = input_params
        self.name = name

        if self.output_dir.name.startswith("boltz_results_"):
            self.output_dir = self.output_dir.rename(
                self.output_dir.parent / f"boltz_{name}"
            )
        self.yaml_input_obj = self.get_input_yaml()
        self.output = self.process_boltz_output()

        self.pae_files = [value["pae"] for value in self.output.values()]
        self.cif_files = [value["cif"] for value in self.output.values()]
        self.pae_to_af3()
        self.af3_pae_files = [value["af3_pae"] for value in self.output.values()]
        self.plddt_files = [value["plddt"] for value in self.output.values()]
        self.pde_files = [value["pde"] for value in self.output.values()]
        self.scores_files = [value["json"] for value in self.output.values()]

    def process_boltz_output(self):
        """
        Boltz 실행 출력을 처리하는 함수입니다.
        """
        file_groups = {}
        for pathway in self.output_dir.rglob("*"):
            number = pathway.stem.split("_model_")[-1]
            if not number.isdigit():
                continue
            number = int(number)

            file_type = pathway.suffix[1:]
            if file_type == FileTypes.NPZ.value:
                file_ = NpzFile(str(pathway))
            elif file_type == FileTypes.CIF.value:
                file_ = CifFile(str(pathway), self.input_params)

            elif file_type == FileTypes.JSON.value:
                file_ = ConfidenceJsonFile(str(pathway))
            else:
                continue
            if number not in file_groups:
                file_groups[number] = [file_]
            else:
                file_groups[number].append(file_)

        model_number_file_type_file = {}
        for model_number, files in file_groups.items():
            intermediate_dict = {}
            for file_ in sorted(files, key=lambda x: x.suffix):
                if file_.pathway.stem.startswith("pae"):
                    intermediate_dict["pae"] = file_
                elif file_.pathway.stem.startswith("plddt"):
                    intermediate_dict["plddt"] = file_
                elif file_.pathway.stem.startswith("pde"):
                    intermediate_dict["pde"] = file_
                elif file_.pathway.suffix == ".cif":
                    file_.name = f"Boltz_{model_number}"
                    file_ = self.update_chain_labels(file_)
                    intermediate_dict["cif"] = file_
                else:
                    intermediate_dict[file_.suffix] = file_

            model_number_file_type_file[model_number] = intermediate_dict

        model_number_file_type_file = {
            key: model_number_file_type_file[key]
            for key in sorted(model_number_file_type_file)
        }
        return model_number_file_type_file

    def add_plddt_to_cif(self):
        """
        PLDDT 점수를 CIF 파일의 B-factor에 추가합니다. Boltz에서는 기본적으로 이 작업을 수행하지 않습니다.

        Returns:
            None

        Raises:
            AssertionError: PLDDT 점수의 길이가 CIF 파일의 잔기 수와 일치하지 않는 경우
        """
        for cif_file, plddt_scores in zip(self.cif_files, self.plddt_files):
            cif_file = self.update_chain_labels(cif_file)
            plddt_score = plddt_scores.data["plddt"]
            if max(plddt_score) <= 1:
                plddt_score = (plddt_score * 100).astype(float)

            chain_lengths = cif_file.chain_lengths(
                mode=ModelCount.RESIDUES, ligand_atoms=True, ptm_atoms=True
            )

            assert sum(chain_lengths.values()) == len(plddt_score), "길이 불일치"

            counter = 0
            for chain in cif_file.model[0]:

                # Boltz는 원자당 리간드 plddt를 수행하므로 별도로 계산해야 합니다.
                ligand = cif_file.check_ligand(chain)

                for residue in chain:
                    for atom in residue:
                        atom.b_iso = plddt_score[counter]
                        if ligand:
                            counter += 1
                    if not ligand:
                        counter += 1

            assert counter == len(plddt_score), "Length mismatch"
            cif_file.update()

    def pae_to_af3(self):
        """
        Boltz의 PAE 데이터를 Alphafold3에서 사용하는 형식으로 변환합니다.

        Returns:
            None
        """
        for i, (pae_file, cif_file) in enumerate(zip(self.pae_files, self.cif_files)):
            pae = Af3Pae.from_boltz(
                pae_file.data,
                cif_file,
            )

            out_name = cif_file.pathway.parent.joinpath(
                cif_file.pathway.stem + "_af3_pae.json"
            )

            pae.to_file(out_name)

            self.output[i]["af3_pae"] = ConfidenceJsonFile(out_name)

    def update_chain_labels(self, cif_file) -> CifFile:
        """
        CIF 파일의 체인 레이블을 업데이트하는 함수입니다.

        Args:
            cif_file (CifFile): 체인 레이블을 업데이트할 CifFile 객체

        """
        cif_file.relabel_chains(
            self.yaml_input_obj.chain_ids, self.yaml_input_obj.id_links
        )
        return cif_file

    def get_input_yaml(self) -> BoltzYaml:
        """
        Boltz 실행에 사용된 입력 yaml 파일을 가져옵니다.

        Returns:
            BoltzYaml: 입력 yaml 파일을 포함하는 객체
        """

        by = BoltzYaml(self.output_dir, create_files=False)
        by.json_to_yaml(self.input_params)

        return by
