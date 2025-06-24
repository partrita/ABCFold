import logging
import subprocess
import sys
from pathlib import Path
from typing import Union

logger = logging.getLogger("logger")


def run_alphafold3(
    input_json: Union[str, Path],
    output_dir: Union[str, Path],
    model_params: Union[str, Path],
    database_dir: Union[str, Path],
    sif_path: Union[str, Path, None],
    interactive: bool = False,
    number_of_models: int = 5,
    num_recycles: int = 10,
) -> bool:
    """
    입력 JSON 파일을 사용하여 Alphafold3를 실행합니다.

    Args:
        input_json (Union[str, Path]): 입력 JSON 파일 경로
        output_dir (Union[str, Path]): 출력 디렉터리 경로
        model_params (Union[str, Path]): 모델 매개변수 경로
        database_dir (Union[str, Path]): 데이터베이스 디렉터리 경로
        sif_path (Union[str, Path, None]): Singularity 이미지 파일 경로
        interactive (bool): True인 경우 Docker 컨테이너를 대화형 모드로 실행합니다.
        number_of_models (int): 생성할 모델 수
        num_recycles (int): 재활용 횟수

    Returns:
        bool: Alphafold3 실행이 성공하면 True, 그렇지 않으면 False

    Raises:
        subprocess.CalledProcessError: Alphafold3 명령이 오류를 반환하는 경우

    """

    input_json = Path(input_json)
    output_dir = Path(output_dir)
    cmd = generate_af3_cmd(
        input_json=input_json,
        output_dir=output_dir,
        model_params=model_params,
        database_dir=database_dir,
        sif_path=sif_path,
        interactive=interactive,
        number_of_models=number_of_models,
        num_recycles=num_recycles,
    )

    logger.info("Running Alphafold3")
    with subprocess.Popen(
        cmd, shell=True, stdout=sys.stdout, stderr=subprocess.PIPE
    ) as p:
        _, stderr = p.communicate()
        if p.returncode != 0:
            logger.error(stderr.decode())
            output_err_file = output_dir / "af3_error.log"
            with open(output_err_file, "w") as f:
                f.write(stderr.decode())
            logger.error("Alphafold3 run failed. Error log is in %s", output_err_file)
            return False

    logger.info("Alphafold3 run complete")
    logger.info("Output files are in %s", output_dir)
    return True


def generate_af3_cmd(
    input_json: Union[str, Path],
    output_dir: Union[str, Path],
    model_params: Union[str, Path],
    database_dir: Union[str, Path],
    sif_path: Union[str, Path, None],
    number_of_models: int = 10,
    num_recycles: int = 5,
    interactive: bool = False,
) -> str:
    """
    Alphafold3 명령을 생성합니다.

    Args:
        input_json (Union[str, Path]): 입력 JSON 파일 경로
        output_dir (Union[str, Path]): 출력 디렉터리 경로
        model_params (Union[str, Path]): 모델 매개변수 경로
        database_dir (Union[str, Path]): 데이터베이스 디렉터리 경로
        sif_path (Union[str, Path, None]): Singularity 이미지 파일 경로
        number_of_models (int): 생성할 모델 수
        num_recycles (int): 재활용 횟수
        interactive (bool): True인 경우 Docker 컨테이너를 대화형 모드로 실행합니다.

    Returns:
        str: Alphafold3 명령
    """
    input_json = Path(input_json)
    output_dir = Path(output_dir)

    if sif_path is not None:
        return f"""
        singularity exec \
        --nv \
        --bind {input_json.parent.resolve()}:/root/af_input \
        --bind {output_dir.resolve()}:/root/af_output \
        --bind {model_params}:/root/models \
        --bind {database_dir}:/root/public_databases \
        {sif_path} \
        python /app/alphafold/run_alphafold.py \
        --json_path=/root/af_input/{input_json.name} \
        --model_dir=/root/models \
        --output_dir=/root/af_output \
        --num_diffusion_samples {number_of_models}\
        --num_recycles {num_recycles}
    """

    else:
        return f"""
        docker run {'-it' if interactive else ''} \
        --volume {input_json.parent.resolve()}:/root/af_input \
        --volume {output_dir.resolve()}:/root/af_output \
        --volume {model_params}:/root/models \
        --volume {database_dir}:/root/public_databases \
        --gpus all \
        alphafold3 \
        python run_alphafold.py \
        --json_path=/root/af_input/{input_json.name} \
        --model_dir=/root/models \
        --output_dir=/root/af_output \
        --num_diffusion_samples {number_of_models}\
        --num_recycles {num_recycles}
        """
