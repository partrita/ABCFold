import argparse
import logging
import sys
from pathlib import Path

logger = logging.getLogger("logger")


def validate_json_file(value):
    """
    입력이 .json 접미사를 가진 JSON 파일인지 확인합니다.
    """
    if not value.endswith(".json"):
        raise argparse.ArgumentTypeError(
            f"입력 파일은 .json 접미사를 가져야 합니다: {value}"
        )
    if not Path(value).exists():
        raise argparse.ArgumentTypeError(f"입력 파일이 존재하지 않습니다: {value}")
    return value


def main_argpase_util(parser):
    parser.add_argument(
        "input_json",
        type=validate_json_file,
        help="AlphaFold3 형식의 입력 JSON 경로",
    )
    parser.add_argument("output_dir", help="출력 디렉터리 경로")
    parser.add_argument(
        "--override",
        help="[선택 사항] 기존 출력 디렉터리가 있는 경우 덮어씁니다.",
        action="store_true",
    )
    parser.add_argument(
        "--output_json",
        help="[선택 사항] 출력 ABCFold JSON 파일의 경로를 지정합니다. \
동일한 입력 기능(예: MSA)으로 ABCFold를 후속 실행하는 데 사용할 수 있습니다.",
    )

    return parser


def mmseqs2_argparse_util(parser):
    parser.add_argument(
        "--mmseqs2",
        action="store_true",
        help="[선택 사항] MSA 생성 및 템플릿 검색에 MMseqs2를 사용합니다 \
(--templates 플래그와 함께 사용하는 경우).",
    )
    parser.add_argument(
        "--mmseqs_database",
        help="[선택 사항] MSA 생성을 위한 데이터베이스 디렉터리입니다. \
로컬 MMseqs2 설치를 사용하는 경우에만 필요합니다.",
    )
    parser.add_argument(
        "--templates", action="store_true", help="[선택 사항] 템플릿 검색 활성화"
    )
    parser.add_argument(
        "--num_templates",
        type=int,
        default=20,
        help="[선택 사항] 사용할 템플릿 수 (기본값: 20)",
    )

    return parser


def custom_template_argpase_util(parser):
    parser.add_argument(
        "--target_id",
        nargs="+",
        help="[조건부 필수] 사용자 정의 템플릿이 관련된 서열의 ID입니다. \
복합체를 모델링하는 경우에만 필요합니다. \
사용자 정의 템플릿 목록을 제공하는 경우 target_id는 \
사용자 정의 템플릿 목록과 길이가 같은 목록이어야 합니다.",
    )
    parser.add_argument(
        "--custom_template",
        nargs="+",
        help="[선택 사항] mmCif 형식의 사용자 정의 템플릿 파일 경로 또는 \
mmCif 형식의 사용자 정의 템플릿 파일 경로 목록입니다. \
사용자 정의 템플릿 목록을 제공하는 경우 사용자 정의 템플릿 체인 목록도 제공해야 합니다.",
    )
    parser.add_argument(
        "--custom_template_chain",
        nargs="+",
        help="[조건부 필수] 사용자 정의 템플릿에서 사용할 체인의 체인 ID입니다. \
다중 체인 템플릿을 사용하는 경우에만 필요합니다. \
사용자 정의 템플릿 목록을 제공하는 경우 사용자 정의 \
템플릿 목록과 길이가 같은 사용자 정의 템플릿 체인 목록을 제공해야 합니다.",
    )

    return parser


def prediction_argparse_util(parser):
    parser.add_argument(
        "--number_of_models",
        type=int,
        default=5,
        help="[선택 사항] 각 방법으로 생성할 모델 수 (기본값: 5)",
    )
    parser.add_argument(
        "--num_recycles",
        type=int,
        default=10,
        help="[선택 사항] 추론 중 사용할 재활용 횟수 (기본값: 10)",
    )
    return parser


def boltz_argparse_util(parser):
    parser.add_argument(
        "-b",
        "--boltz",
        action="store_true",
        help="Boltz 실행",
    )
    if "--save_input" not in parser._option_string_actions:
        parser.add_argument(
            "--save_input",
            action="store_true",
            help="입력 JSON 파일 저장",
            default=False,
        )

    return parser


def chai_argparse_util(parser):
    parser.add_argument(
        "-c",
        "--chai1",
        action="store_true",
        help="Chai-1 실행",
    )
    return parser


def alphafold_argparse_util(parser):
    parser.add_argument(
        "--database",
        help="[선택 사항] MSA 생성을 위한 데이터베이스 디렉터리입니다. \
내장 AlphaFold3 MSA 생성을 사용하는 경우에만 필요합니다.",
        dest="database_dir",
        default=None,
    )

    parser.add_argument(
        "--model_params",
        help="[필수] AlphaFold3 모델 매개변수를 포함하는 디렉터리",
        default=None,
    )

    parser.add_argument(
        "--sif_path",
        help="[조건부 필수] Singularity를 사용하는 경우 AlphaFold3의 sif 이미지 경로",
        default=None,
    )

    parser.add_argument(
        "-a",
        "--alphafold3",
        action="store_true",
        help="Alphafold3 실행",
    )

    parser.add_argument(
        "--use_af3_template_search",
        action="store_true",
        help="자체 사용자 정의 MSA를 제공하거나 `--mmseqs2`를 실행한 경우, \
Alphafold3가 템플릿을 검색하도록 허용합니다.",
    )

    return parser


def visuals_argparse_util(parser):
    parser.add_argument(
        "--no_visuals",
        action="store_true",
        help="[선택 사항] 출력 페이지를 생성하지 않습니다. 디스플레이가 없는 \
클러스터에서 실행하는 데 가장 적합합니다.",
    )

    parser.add_argument(
        "--no_server",
        action="store_true",
        help="[선택 사항] 결과를 보기 위해 로컬 서버를 시작하지 않습니다. 출력 \
페이지는 여전히 생성되며 출력 디렉터리에서 액세스할 수 있습니다.",
    )
    return parser


def raise_argument_errors(args):
    if not args.alphafold3 and not args.boltz and not args.chai1:
        logger.info(
            "AlphaFold3, Boltz 또는 Chai-1이 선택되지 않았습니다. 기본적으로 AlphaFold3을 \
실행합니다."
        )
        args.alphafold3 = True

    if (
        args.alphafold3
        and (not args.model_params or not Path(args.model_params).exists())
        and not args.mmseqs2
    ):
        logger.error(f"모델 매개변수 디렉터리를 찾을 수 없습니다: {args.model_params}")
        sys.exit(1)

    if args.templates and not args.mmseqs2 and not args.alphafold3:
        logger.error("--templates 플래그는 MMseqs2 또는 Alphafold3 없이 사용할 수 없습니다.")
        sys.exit(1)

    if (
        args.templates
        and args.alphafold3
        and not args.mmseqs2
        and not args.use_af3_template_search
    ):
        # --templates가 설정된 경우 Alphafold3와 함께 템플릿이 사용되도록 보장
        args.use_af3_template_search = True

    if args.custom_template_chain and not args.custom_template:
        logger.error("사용자 정의 템플릿 없이 사용자 정의 템플릿 체인이 제공되었습니다.")
        sys.exit(1)

    if args.use_af3_template_search and not args.alphafold3:
        logger.error(
            "Alphafold3를 실행하지 않고는 Alphafold3 템플릿 검색을 사용할 수 없습니다."
        )
        sys.exit(1)

    if args.num_templates < 1:
        logger.error("템플릿 수는 0보다 커야 합니다.")
        sys.exit(1)

    if args.num_recycles < 1:
        logger.error("재활용 횟수는 0보다 커야 합니다.")
        sys.exit(1)

    if args.number_of_models < 1:
        logger.error("모델 수는 0보다 커야 합니다.")
        sys.exit(1)

    return args
