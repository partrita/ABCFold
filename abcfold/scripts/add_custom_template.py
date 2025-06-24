#!/usr/bin/env python

import json
import logging
import os

from abcfold.argparse_utils import custom_template_argpase_util
from abcfold.scripts.abc_script_utils import get_custom_template

logger = logging.getLogger("logger")


def run_custom_template(
    input_json,
    target_id,
    custom_template,
    custom_template_chain,
    output_json=None,
    to_file=True,
):
    """
    AlphaFold3 입력 JSON 객체에 사용자 정의 템플릿을 추가합니다.

    Args:
        input_json (str or Path): 입력 AlphaFold3 JSON 파일 경로.
        target_id (list or None): 사용자 정의 템플릿에 대한 대상 ID 목록.
        custom_template (list): 사용자 정의 템플릿 CIF 파일 경로 목록.
        custom_template_chain (list): 사용자 정의 템플릿에서 사용할 체인 ID 목록.
        output_json (str or Path, optional): 수정된 JSON을 저장할 경로.
                                            None이면 input_json을 덮어씁니다. 기본값은 None입니다.
        to_file (bool, optional): 출력을 파일에 저장할지 여부. 기본값은 True입니다.

    Returns:
        dict: 수정된 AlphaFold3 JSON 객체.

    Raises:
        FileNotFoundError: 사용자 정의 템플릿 파일을 찾을 수 없는 경우.
        ValueError: 제공된 인수에 불일치가 있는 경우
                    (예: 목록 길이 불일치, 다중 서열에 대한 target_id 누락).
    """
    af3_json = json.load(open(input_json))

    for sequence in af3_json["sequences"]:
        if "protein" not in sequence:
            continue

        for template in custom_template:
            if not os.path.exists(template):
                msg = f"사용자 정의 템플릿 파일 {template}을(를) 찾을 수 없습니다."
                logger.critical(msg)
                raise FileNotFoundError()
            # 단백질 서열에만 템플릿을 추가할 수 있으므로 입력 JSON에
            # 여러 단백질 서열이 있는지 확인합니다.
            if (
                len([x for x in af3_json["sequences"] if "protein" in x.keys()]) > 1
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
            custom_templates = zip(target_id, custom_template, custom_template_chain)
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
            custom_templates = zip(target_ids, custom_template, custom_template_chain)

        for i in custom_templates:
            tid, c_tem, c_tem_chn = i
            sequence = get_custom_template(
                sequence,
                tid,
                c_tem,
                c_tem_chn,
            )

    if to_file:
        if not output_json:
            output_json = input_json

        with open(output_json, "w") as f:
            json.dump(af3_json, f)

    return af3_json


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="AlphaFold 입력 JSON에 사용자 정의 템플릿 추가"
    )

    parser.add_argument("--input_json", help="입력 AlphaFold3 JSON 파일")
    parser.add_argument("--output_json", help="Output alphafold3 json file")
    parser = custom_template_argpase_util(parser)

    args = parser.parse_args()

    run_custom_template(  # pragma: no cover
        args.input_json,
        args.target_id,
        args.custom_template,
        args.custom_template_chain,
        output_json=args.output_json,
        to_file=True,
    )


if __name__ == "__main__":  # pragma: no cover
    main()
