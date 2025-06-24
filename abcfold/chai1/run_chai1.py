import logging
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Union

from abcfold.chai1.af3_to_chai import ChaiFasta
from abcfold.chai1.check_install import check_chai1

logger = logging.getLogger("logger")


def run_chai(
    input_json: Union[str, Path],
    output_dir: Union[str, Path],
    save_input: bool = False,
    test: bool = False,
    number_of_models: int = 5,
    num_recycles: int = 10,
    use_templates_server: bool = False,
    template_hits_path: Path | None = None,
) -> bool:
    """
    입력 JSON 파일을 사용하여 Chai-1을 실행합니다.

    Args:
        input_json (Union[str, Path]): 입력 JSON 파일 경로
        output_dir (Union[str, Path]): 출력 디렉터리 경로
        save_input (bool): True인 경우 입력 fasta 파일과 MSA를 출력 디렉터리에 저장합니다.
        test (bool): True인 경우 테스트 명령을 실행합니다.
        number_of_models (int): 생성할 모델 수
        num_recycles (int): 트렁크 재활용 횟수
        use_templates_server (bool): True인 경우 서버의 템플릿을 사용합니다.
        template_hits_path (Path | None): 템플릿 히트 m8 파일 경로

    Returns:
        bool: Chai-1 실행이 성공하면 True, 그렇지 않으면 False

    """
    input_json = Path(input_json)
    output_dir = Path(output_dir)

    # Chai-1 설치 여부 확인
    logger.debug("Checking if Chai-1 is installed")
    check_chai1()

    with tempfile.TemporaryDirectory() as temp_dir:
        working_dir = Path(temp_dir)
        chai_output_dir = output_dir
        if save_input:
            logger.info("Saving input fasta file and msa to the output directory")
            working_dir = output_dir
            working_dir.mkdir(parents=True, exist_ok=True)
            chai_output_dir = output_dir / "chai_output"

        chai_fasta = ChaiFasta(working_dir)
        chai_fasta.json_to_fasta(input_json)

        out_fasta = chai_fasta.fasta
        msa_dir = chai_fasta.working_dir
        out_constraints = chai_fasta.constraints

        cmd = (
            generate_chai_command(
                out_fasta,
                msa_dir,
                out_constraints,
                chai_output_dir,
                number_of_models,
                num_recycles=num_recycles,
                use_templates_server=use_templates_server,
                template_hits_path=template_hits_path,
            )
            if not test
            else generate_chai_test_command()
        )

        logger.info("Running Chai-1")
        with subprocess.Popen(
            cmd,
            stdout=sys.stdout,
            stderr=subprocess.PIPE,
        ) as proc:
            _, stderr = proc.communicate()
            if proc.returncode != 0:
                if proc.stderr:
                    if chai_output_dir.exists():
                        output_err_file = chai_output_dir / "chai_error.log"
                    else:
                        output_err_file = chai_output_dir.parent / "chai_error.log"
                    with open(output_err_file, "w") as f:
                        f.write(stderr.decode())
                    logger.error(
                        "Chai-1 run failed. Error log is in %s", output_err_file
                    )
                else:
                    logger.error("Chai-1 run failed")
                return False

        logger.info("Chai-1 run complete")
        return True


def generate_chai_command(
    input_fasta: Union[str, Path],
    msa_dir: Union[str, Path],
    input_constraints: Union[str, Path],
    output_dir: Union[str, Path],
    number_of_models: int = 5,
    num_recycles: int = 10,
    use_templates_server: bool = False,
    template_hits_path: Path | None = None,
) -> list:
    """
    Chai-1 명령을 생성합니다.

    Args:
        input_fasta (Union[str, Path]): 입력 fasta 파일 경로
        msa_dir (Union[str, Path]): MSA 디렉터리 경로
        input_constraints (Union[str, Path]): 입력 제약 조건 파일 경로
        output_dir (Union[str, Path]): 출력 디렉터리 경로
        number_of_models (int): 생성할 모델 수
        num_recycles (int): 트렁크 재활용 횟수
        use_templates_server (bool): True인 경우 서버의 템플릿을 사용합니다.
        template_hits_path (Path | None): 템플릿 히트 m8 파일 경로

    Returns:
        list: Chai-1 명령

    """

    chai_exe = Path(__file__).parent / "chai.py"
    cmd = ["python", str(chai_exe), "fold", str(input_fasta)]

    if Path(msa_dir).exists():
        cmd += ["--msa-directory", str(msa_dir)]
    if Path(input_constraints).exists():
        cmd += ["--constraint-path", str(input_constraints)]

    cmd += ["--num-diffn-samples", str(number_of_models)]
    cmd += ["--num-trunk-recycles", str(num_recycles)]

    assert not (
        use_templates_server and template_hits_path
    ), "Cannot specify both templates server and path"

    if shutil.which("kalign") is None and (use_templates_server or template_hits_path):
        logger.warning(
            "kalign not found, skipping template search kalign is required. \
Please install kalign to use templates with Chai-1."
        )
    else:
        if use_templates_server:
            cmd += ["--use-templates-server"]
        if template_hits_path:
            cmd += ["--template-hits-path", str(template_hits_path)]

    cmd += [str(output_dir)]

    return cmd


def generate_chai_test_command() -> list:
    """
    Chai-1 테스트 명령을 생성합니다.

    Args:
        None

    Returns:
        list: Chai-1 테스트 명령
    """
    chai_exe = Path(__file__).parent / "chai.py"
    return [
        "python",
        chai_exe,
        "fold",
        "--help",
    ]
