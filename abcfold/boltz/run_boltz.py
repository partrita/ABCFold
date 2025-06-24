import logging
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Union

from abcfold.boltz.af3_to_boltz import BoltzYaml
from abcfold.boltz.check_install import check_boltz

logger = logging.getLogger("logger")


def run_boltz(
    input_json: Union[str, Path],
    output_dir: Union[str, Path],
    save_input: bool = False,
    test: bool = False,
    number_of_models: int = 5,
    num_recycles: int = 10,
) -> bool:
    """
    입력 JSON 파일을 사용하여 Boltz를 실행합니다.

    Args:
        input_json (Union[str, Path]): 입력 JSON 파일 경로
        output_dir (Union[str, Path]): 출력 디렉터리 경로
        save_input (bool): True인 경우 입력 yaml 파일과 MSA를 출력 디렉터리에 저장합니다.
        test (bool): True인 경우 테스트 명령을 실행합니다.
        number_of_models (int): 생성할 모델 수
        num_recycles (int): 재활용 횟수

    Returns:
        bool: Boltz 실행이 성공하면 True, 그렇지 않으면 False

    Raises:
        subprocess.CalledProcessError: Boltz 명령이 오류를 반환하는 경우

    """
    input_json = Path(input_json)
    output_dir = Path(output_dir)

    # Boltz 설치 여부 확인
    logger.debug("Checking if boltz is installed")
    check_boltz()

    with tempfile.TemporaryDirectory() as temp_dir:
        working_dir = Path(temp_dir)
        if save_input:
            logger.info("Saving input yaml file and msa to the output directory")
            working_dir = output_dir

        boltz_yaml = BoltzYaml(working_dir)
        boltz_yaml.json_to_yaml(input_json)
        out_file = working_dir.joinpath(f"{input_json.stem}.yaml")

        boltz_yaml.write_yaml(out_file)
        logger.info("Running Boltz")
        cmd = (
            generate_boltz_command(out_file, output_dir, number_of_models, num_recycles)
            if not test
            else generate_boltz_test_command()
        )

        with subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ) as proc:
            stdout = ""
            if proc.stdout:
                for line in proc.stdout:
                    sys.stdout.write(line.decode())
                    sys.stdout.flush()
                    stdout += line.decode()
            _, stderr = proc.communicate()
            if proc.returncode != 0:
                if proc.stderr:
                    logger.error(stderr.decode())
                    output_err_file = output_dir / "boltz_error.log"
                    with open(output_err_file, "w") as f:
                        f.write(stderr.decode())
                    logger.error(
                        "Boltz run failed. Error log is in %s", output_err_file
                    )
                else:
                    logger.error("Boltz run failed")
                return False
            elif "WARNING: ran out of memory" in stdout:
                logger.error("Boltz ran out of memory")
                return False

        logger.info("Boltz run complete")
        logger.info("Output files are in %s", output_dir)
        return True


def generate_boltz_command(
    input_yaml: Union[str, Path],
    output_dir: Union[str, Path],
    number_of_models: int = 5,
    num_recycles: int = 10,
) -> list:
    """
    Boltz 명령을 생성합니다.

    Args:
        input_yaml (Union[str, Path]): 입력 YAML 파일 경로
        output_dir (Union[str, Path]): 출력 디렉터리 경로
        number_of_models (int): 생성할 모델 수
        num_recycles (int): 재활용 횟수

    Returns:
        list: Boltz 명령
    """
    return [
        "boltz",
        "predict",
        str(input_yaml),
        "--out_dir",
        str(output_dir),
        "--override",
        "--write_full_pae",
        "--write_full_pde",
        "--diffusion_samples",
        str(number_of_models),
        "--recycling_steps",
        str(num_recycles),
    ]


def generate_boltz_test_command() -> list:
    """
    Boltz에 대한 테스트 명령을 생성합니다.

    Args:
        None

    Returns:
        list: Boltz 테스트 명령
    """

    return [
        "boltz",
        "predict",
        "--help",
    ]
