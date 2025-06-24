import logging
from pathlib import Path
from typing import Union

from abcfold.chai1.af3_to_chai import ChaiFasta
from abcfold.output.file_handlers import (CifFile, ConfidenceJsonFile,
                                          FileTypes, NpyFile, NpzFile)
from abcfold.output.utils import Af3Pae

logger = logging.getLogger("logger")


class ChaiOutput:
    def __init__(
        self,
        chai_output_dir: Union[str, Path],
        input_params: dict,
        name: str,
        save_input: bool = False,
    ):
        """
        Chai-1 실행 출력을 처리하는 객체입니다.

        Args:
            chai_output_dir (Union[str, Path]): Chai-1 출력 디렉터리 경로
            input_params (dict): Chai-1 실행에 사용된 입력 매개변수를 포함하는 딕셔너리
            name (str): Chai-1 실행에 지정된 이름
            save_input (bool): True인 경우 Chai-1이 save_input 플래그로 실행되었습니다.

        Attributes:
            input_params (dict): Chai-1 실행에 사용된 입력 매개변수를 포함하는 딕셔너리
            output_dir (Path): Chai-1 출력 디렉터리 경로
            name (str): Chai-1 실행에 지정된 이름
            output (dict): Chai-1 출력 디렉터리의 내용을 처리한 출력을 포함하는 딕셔너리입니다.
                           딕셔너리 구조는 다음과 같습니다:

            {
                1: {
                    "pae": NpzFile,
                    "cif": CifFile,
                    "scores": NpyFile
                },
                2: {
                    "pae": NpzFile,
                    "cif": CifFile,
                    "scores": NpyFile
                },
                ...
            }
            pae_files (list): PAE 데이터를 포함하는 정렬된 NpzFile 객체 목록
            cif_files (list): CIF 데이터를 포함하는 정렬된 CifFile 객체 목록
            scores_files (list): 점수 데이터를 포함하는 정렬된 NpyFile 객체 목록

        """
        self.input_params = input_params
        self.output_dir = Path(chai_output_dir)
        self.name = name
        self.save_input = save_input

        if not self.output_dir.name.startswith("chai1_" + self.name):
            self.output_dir = self.output_dir.rename(
                self.output_dir.parent / f"chai1_{self.name}"
            )

        self.input_fasta = self.get_input_fasta()

        self.output = self.process_chai_output()
        self.pae_files = [
            value["pae"] for value in self.output.values() if "pae" in value
        ]
        self.cif_files = [
            value["cif"] for value in self.output.values() if "cif" in value
        ]
        self.pae_to_af3()
        self.scores_files = [
            value["scores"] for value in self.output.values() if "scores" in value
        ]
        self.af3_pae_files = [
            value["af3_pae"] for value in self.output.values() if "af3_pae" in value
        ]

    def process_chai_output(self):
        file_groups = {}

        if self.save_input:
            self.output_dir = self.output_dir / "chai_output"

        for pathway in self.output_dir.iterdir():
            number = pathway.stem.split("model_idx_")[-1]
            if number.isdigit():
                number = int(number)

            file_type = pathway.suffix[1:]

            if file_type == FileTypes.NPZ.value:
                file_ = NpzFile(str(pathway))

            elif file_type == FileTypes.CIF.value:
                file_ = CifFile(str(pathway), self.input_params)
                file_ = self.update_chain_labels(file_)

            elif file_type == FileTypes.NPY.value:
                file_ = NpyFile(str(pathway))
            else:
                continue

            if isinstance(number, str):
                number = -1

            if number not in file_groups:
                file_groups[number] = [file_]
            else:
                file_groups[number].append(file_)

        model_number_file_type_file = {}
        for model_number, files in file_groups.items():
            intermediate_dict = {}
            for file_ in sorted(files, key=lambda x: x.suffix):
                if file_.pathway.stem.startswith("scores.model"):
                    intermediate_dict["scores"] = file_
                elif file_.pathway.stem.startswith("pred.model"):
                    file_.name = f"Chai-1_{model_number}"
                    # Chai cif는 pae-viewer에서 인식되지 않으므로 로드 후 저장합니다.
                    file_.to_file(file_.pathway)
                    intermediate_dict["cif"] = file_
                elif file_.pathway.stem.startswith("pae_scores"):
                    intermediate_dict["pae"] = file_

            model_number_file_type_file[model_number] = intermediate_dict

        model_number_file_type_file = {
            model_number: model_number_file_type_file[model_number]
            for model_number in sorted(model_number_file_type_file)
        }

        return model_number_file_type_file

    def pae_to_af3(self) -> None:
        """
        Chai-1 PAE 데이터를 AlphaFold3에서 예상하는 형식으로 변환합니다.

        """

        pae_file = self.pae_files[-1]
        for i, cif_file in enumerate(self.cif_files):
            pae = Af3Pae.from_chai1(
                pae_file.data[i],
                cif_file,
            )

            out_name = self.output_dir.joinpath(cif_file.pathway.stem + "_af3_pae.json")
            pae.to_file(out_name)

            self.output[i]["af3_pae"] = ConfidenceJsonFile(out_name)

    def get_input_fasta(self) -> ChaiFasta:
        """
        Chai-1 실행에 사용된 입력 fasta 파일을 가져오는 함수입니다.

        Returns:
            ChaiFasta: 입력 fasta 파일을 포함하는 ChaiFasta 객체

        """

        ch = ChaiFasta(self.output_dir, create_files=False)
        ch.json_to_fasta(self.input_params)

        return ch

    def update_chain_labels(self, cif_file: CifFile) -> CifFile:
        """
        CIF 파일의 체인 레이블을 업데이트하는 함수입니다.

        Args:
            cif_file (CifFile): 체인 레이블을 업데이트할 CifFile 객체

        """

        cif_file.relabel_chains(self.input_fasta.chain_ids)
        cif_file.to_file(cif_file.pathway)
        return CifFile(cif_file.pathway, self.input_params)
