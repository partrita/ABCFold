import http.server
import textwrap
from itertools import groupby
from operator import itemgetter
from pathlib import Path
from typing import Dict, Union

import numpy as np
from Bio.SeqUtils import seq1
from jinja2 import Environment, FileSystemLoader

from abcfold.output.alphafold3 import AlphafoldOutput
from abcfold.output.boltz import BoltzOutput
from abcfold.output.chai import ChaiOutput
from abcfold.output.file_handlers import ConfidenceJsonFile, NpzFile
from abcfold.plots.pae_plot import create_pae_plots
from abcfold.plots.plddt_plot import plot_plddt

PORT = 8000


def get_plddt_regions(plddts: Union[np.ndarray, list]) -> dict:
    """
    모델의 pLDDT 영역을 가져옵니다.
    """
    if not isinstance(plddts, np.ndarray):
        plddts = np.array(plddts)

    regions = {}
    # None 값을 -1로 바꿉니다.
    plddts = np.where(plddts is None, -1, plddts)

    v_low = np.where((0 <= plddts) & (plddts <= 50))[0]
    regions["v_low"] = get_regions_helper(v_low)
    low = np.where((plddts > 50) & (plddts < 70))[0]
    regions["low"] = get_regions_helper(low)
    confident = np.where((plddts >= 70) & (plddts < 90))[0]
    regions["confident"] = get_regions_helper(confident)
    v_confident = np.where(plddts >= 90)[0]
    regions["v_high"] = get_regions_helper(v_confident)

    return regions


def get_regions_helper(indices):
    """
    인덱스에서 영역을 가져옵니다.
    """
    regions = []
    for _, g in groupby(enumerate(indices), lambda x: x[0] - x[1]):
        group = map(itemgetter(1), g)
        group = list(map(int, group))
        regions.append((group[0], group[-1]))
    return regions


def get_model_sequence_data(cif_objs) -> dict:
    """
    모델의 각 체인 및 리간드에 대한 서열을 가져옵니다. 플로팅을 위해 내부적으로 사용됩니다.

    Args:
        cif_objs : CifFile 객체 목록

    Returns:
        dict : 체인 ID 및 서열 데이터
    """
    sequence_data: dict = {}
    for cif_obj in cif_objs:
        sequence_data_ = {}
        for chain in cif_obj.get_chains():
            if cif_obj.check_ligand(chain):
                if chain.id not in sequence_data_:
                    sequence_data_[chain.id] = ""
                sequence_data_[chain.id] += "".join(
                    [atom.id[0] for residue in chain for atom in residue]
                )
            elif cif_obj.check_other(chain, ["dna"]):
                sequence_data_[chain.id] = "".join(
                    [residue.get_resname()[-1] for residue in chain]
                )
            elif cif_obj.check_other(chain, ["rna"]):
                sequence_data_[chain.id] = "".join(
                    [residue.get_resname() for residue in chain]
                )
            else:

                sequence_data_[chain.id] = "".join(
                    [seq1(residue.get_resname()) for residue in chain]
                )
        sequence_data = {
            chain_id: sorted(
                [sequence_data_[chain_id], sequence_data.get(chain_id, "")],
                reverse=True,
            )[0]
            for chain_id in sequence_data_
        }

    return sequence_data


def get_model_data(model, plot_dict, method, plddt_scores, score_file, output_dir):
    """
    출력 페이지에 대한 모델 데이터를 가져옵니다.

    Args:
        model (CifFile): 모델 객체
        plot_dict (dict): 플롯 딕셔너리
        method (str): 모델 생성에 사용된 방법
        plddt_scores : pLDDT 점수
        score_file (str): 모델 점수를 포함하는 파일 경로
        output_dir (Path): 출력 디렉터리 경로
    """
    regions = get_plddt_regions(plddt_scores)
    ptm_score, iptm_score = parse_scores(score_file)
    model_path = Path(model.pathway).relative_to(output_dir)
    model_data = {
        "model_id": model.name,
        "model_source": method,
        "model_path": model_path.as_posix(),
        "plddt_regions": regions,
        "avg_plddt": model.average_plddt,
        "h_score": model.h_score,
        "ptm_score": ptm_score,
        "iptm_score": iptm_score,
        "residue_clashes": model.clashes_residues,
        "atom_clashes": model.clashes,
        "pae_path": Path(plot_dict[model.pathway.as_posix()])
        .relative_to(output_dir)
        .as_posix(),
    }
    return model_data


class NoCacheHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header(
            "Cache-Control", "no-store, no-cache, must-revalidate, max-age=0"
        )
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        super().end_headers()


def plots(outputs: list, output_dir: Path):
    """
    다른 프로그램의 출력에 대한 플롯을 생성합니다.

    Args:
        outputs (list): 출력 객체 목록
        output_dir (Path): 출력 디렉터리 경로

    """
    pathway_plots = create_pae_plots(outputs, output_dir=output_dir)
    plddt_plot_input: Dict[str, list] = get_all_cif_files(outputs)

    plot_plddt(plddt_plot_input, output_name=output_dir.joinpath("plddt_plot.html"))

    pathway_plots["plddt"] = str(output_dir.joinpath("plddt_plot.html").resolve())

    return pathway_plots


def render_template(in_file_path, out_file_path, **kwargs):
    """
    주어진 파일을 키워드 인수로 템플릿화합니다.

    Args:
        in_file_path (Path): 템플릿 경로.
        out_file_path (Path): 템플릿화된 파일을 출력할 경로.
        **kwargs (dict): 템플릿팅에 사용할 변수.
    """
    env = Environment(
        loader=FileSystemLoader(in_file_path.parent), keep_trailing_newline=True
    )
    template = env.get_template(in_file_path.name)
    # kwargs는 템플릿에서 변수로 나타남
    output = template.render(**kwargs)
    with open(str(out_file_path), "w") as f:
        f.write(output)


def output_open_html_script(file_out: str, port: int = 8000):
    """
    기본 웹 브라우저에서 출력 HTML 파일을 여는 파이썬 스크립트를 만듭니다.

    Args:
        file_out (str): 출력 스크립트 경로
        port (int): 서버를 실행할 포트
    """

    script = f"""
    import http.server
    import socketserver
    import webbrowser
    import sys

    PORT = {port}

    class NoCacheHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
        def end_headers(self):
            self.send_header("Cache-Control",
                            "no-store, no-cache, must-revalidate, max-age=0")
            self.send_header("Pragma", "no-cache")
            self.send_header("Expires", "0")
            super().end_headers()

    try:
        with socketserver.TCPServer(("", PORT),
                                    NoCacheHTTPRequestHandler) as httpd:
            print(
                f"Serving at port {PORT}: http://localhost:{PORT}/index.html"
                )
            print("Press Ctrl+C to stop the server")
            webbrowser.open(f"http://localhost:{PORT}/index.html")
            httpd.serve_forever()
    except KeyboardInterrupt:
        print("Server stopped")
        httpd.server_close()
        sys.exit(0)
    """

    script = textwrap.dedent(script)
    with open(file_out, "w") as f:
        f.write(script)


def get_all_cif_files(outputs) -> Dict[str, list]:
    """
    모든 CIF 파일을 가져옵니다.

    Args:
        outputs (list): 출력 객체 목록

    Returns:
        Dict[str, list]: 메서드 이름과 CIF 파일 객체 목록을 매핑하는 딕셔너리
    """
    method_cif_objs: Dict[str, list] = {}

    for output in outputs:
        if isinstance(output, AlphafoldOutput):
            for seed in output.seeds:
                if "Alphafold3" not in method_cif_objs:
                    method_cif_objs["Alphafold3"] = []
                method_cif_objs["Alphafold3"].extend(output.cif_files[seed])
        elif isinstance(output, BoltzOutput):

            method_cif_objs["Boltz"] = output.cif_files
        elif isinstance(output, ChaiOutput):
            method_cif_objs["Chai-1"] = output.cif_files

    return method_cif_objs


def parse_scores(score_file: Union[ConfidenceJsonFile, NpzFile]) -> tuple:
    """
    점수 파일에서 점수를 파싱합니다.

    Args:
        score_file (Union[ConfidenceJsonFile, NpzFile]): 점수 파일 객체.

    Returns:
        tuple: ptm_score와 iptm_score를 float으로 포함하는 튜플.
    """
    ptm_score = None
    iptm_score = None

    if isinstance(score_file, ConfidenceJsonFile):
        data = score_file.load_json_file()
        if "ptm" in data and "iptm" in data:
            ptm_score = round(float(data["ptm"]), 2)
            iptm_score = round(float(data["iptm"]), 2)
    elif isinstance(score_file, NpzFile):
        data = score_file.load_npz_file()
        if "ptm" in data and "iptm" in data:
            ptm_score = round(float(data["ptm"]), 2)
            iptm_score = round(float(data["iptm"]), 2)

    return ptm_score, iptm_score
