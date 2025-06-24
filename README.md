# ABCFold

![빌드 상태](https://github.com/rigdenlab/ABCFold/actions/workflows/python-package.yml/badge.svg)
![커버리지](https://raw.githubusercontent.com/rigdenlab/ABCFold/refs/heads/main/.blob/coverage.svg)

MMseqs2 다중 서열 정렬(MSA) 및 사용자 정의 템플릿을 사용하여 AlphaFold3, Boltz 및 Chai-1을 실행하는 스크립트입니다.

## 목차
- [설치](#설치)
- [사용법](#사용법)
- [일반적인 문제](#일반적인-문제)
- [기여](#기여)

## 설치

가상 환경 또는 conda / micromamba 환경에 이 패키지를 설치하는 것을 권장합니다. Python 3.11이 권장되지만, Python 3.9 이상에서도 작동해야 합니다.

conda/micromamba 환경을 설정하려면 다음을 실행하십시오:
```bash
conda create -n abcfold python=3.11
conda activate abcfold
```

또는

```bash
micromamba env create -n abcfold python=3.11
micromamba activate abcfold
```

PyPI에서 패키지를 설치하려면 다음을 실행하십시오:

```bash
python -m pip install abcfold
```

또는 소스에서 패키지를 설치하려면 먼저 저장소를 복제한 다음 실행하십시오:

```bash
python -m pip install .
```

## 개발

이 패키지 개발에 참여하고 싶다면 다음을 실행하여 개발 종속성을 설치할 수 있습니다:

```bash
python -m pip install -e .
python -m pip install -r requirements-dev.txt
python -m pre_commit install
```

## 사용법

### ABCfold 실행

ABCFold는 Alphafold3, Boltz 및 Chai-1을 순차적으로 실행합니다. 이 프로그램은 Alphafold3 형식의 JSON 입력을 받습니다 (이 형식 지정 방법에 대한 전체 지침은 [여기](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md)를 클릭하십시오). 예제 JSON은 다음과 같습니다:

```json
{
  "name": "2PV7",
  "sequences": [
    {
      "protein": {
        "id": ["A", "B"],
        "sequence": "GMRESYANENQFGFKTINSDIHKIVIVGGYGKLGGLFARYLRASGYPISILDREDWAVAESILANADVVIVSVPINLTLETIERLKPYLTENMLLADLTSVKREPLAKMLEVHTGAVLGLHPMFGADIASMAKQVVVRCDGRFPERYEWLLEQIQIWGAKIYQTNATEHDHNMTYIQALRHFSTFANGLHLSKQPINLANLLALSSPIYRLELAMIGRLFAQDAELYADIIMDKSENLAVIETLKQTYDEALTFFENNDRQGFIDAFHKVRDWFGDYSEQFLKESRQLLQQANDLKQG"
      }
    }
  ],
  "modelSeeds": [1],
  "dialect": "alphafold3",
  "version": 1
}
```

시스템에 AlphaFold3가 설치되어 있고 ([여기](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md) 지침) 모델 매개변수를 확보했는지 확인하십시오. Boltz 및 Chai-1은 런타임 시 설치됩니다.

대부분의 작업에서 ABCFold는 다음과 같이 실행할 수 있습니다:
```bash
abcfold <input_json>  <output_dir> -abc --mmseqs2 --model_params <path_to_af3_model_params>
```
> [!NOTE]
> `--model_params`는 첫 실행 후 저장되므로, 이후 ABCFold 작업에는 이 플래그가 필요하지 않습니다.

> [!NOTE]
> AlphaFold3 JACKHMMER MSA 검색으로 ABCFold를 실행하려면 `--mmseqs2` 플래그를 제거하고 AlphaFold3 데이터베이스가 포함된 디렉터리 경로와 함께 `--database` 플래그를 제공해야 합니다.
> `--database` 경로도 첫 실행 후 저장되며 이후 ABCFold 작업에는 필요하지 않습니다.

>[!WARNING]
> `--mmseqs2` 플래그를 사용하면 AlphaFold3가 pairedMSA 정보 없이 실행됩니다. 이것이 대상에 중요한 경우(예: 복합체 모델링), pairedMSA가 자동으로 생성되므로 AlphaFold3 JACKHMMER MSA 검색을 실행하는 것이 좋습니다.

>[!WARNING]
> 새로 설치하는 경우 `--model_params`와 `--database`를 다시 제공해야 합니다.

그러나 템플릿 사용이나 생성할 모델 수와 같은 런타임 옵션을 추가하기 위해 다음 플래그를 사용할 수 있습니다.

#### 주요 인수
- `<input_json>`: 입력 AlphaFold3 JSON 파일 경로.
- `<output_dir>`: 출력 디렉터리 경로.
- `-a`, `-b`, `-c` (`--alphafold3`, `--boltz`,`--chai1`): 각각 Alphafold3, Boltz, Chai-1을 실행하는 플래그. 이 플래그 중 아무것도 제공되지 않으면 기본적으로 Alphafold3가 실행됩니다.
- `--mmseqs2`: [선택 사항] MMseqs2 MSA 및 템플릿(지정된 경우)을 사용하는 플래그.
- `--mmseqs_database`: [선택 사항] 로컬 MMSeqs2 사본에서 사용하는 데이터베이스 경로. mmseqs가 설치되어 있다면 이 플래그를 포함하면 MMseqs2를 로컬에서 실행할 수 있습니다.
- `--override`: [선택 사항] 기존 출력 디렉터리를 덮어쓰는 플래그.
- `--save_input`: [선택 사항] 입력 JSON 파일을 출력 디렉터리에 저장하는 플래그.

#### Alphafold3 인수

- `--model_params`: AlphaFold3 모델 매개변수가 포함된 디렉터리 경로.
- `--database`: [선택 사항] AlphaFold3 데이터베이스가 포함된 디렉터리 경로 #참고: `--mmseqs2` 플래그를 사용하는 경우에는 사용되지 않습니다.
- `--sif_path`: [선택 사항] Docker 대신 AlphaFold3 특이점을 사용하는 경우 sif 파일 경로
- `--use_af3_template_search`[선택 사항] 사용자 정의 MSA를 제공하거나 `--mmseqs2`를 실행한 경우 Alphafold3가 템플릿을 검색하도록 허용합니다.

#### 템플릿 인수

- `--templates`: 템플릿 검색을 활성화하는 플래그
- `--num_templates`: [선택 사항] 사용할 템플릿 수 (기본값: 20)

- `--custom_template`: [선택 사항] mmCIF 형식의 사용자 정의 템플릿 파일 경로 또는 사용자 정의 템플릿 목록. 사용자 정의 템플릿 인수 사용 방법에 대한 자세한 설명은 아래 시각화 인수에서 찾을 수 있습니다.
- `--custom_template_chain`: [조건부 필수] 사용자 정의 템플릿에서 사용할 체인의 체인 ID. 다중 체인 템플릿을 사용하는 경우에만 필요합니다. 사용자 정의 템플릿 목록을 제공하는 경우 사용자 정의 템플릿 체인 목록을 제공해야 합니다.
- `--target_id`: [조건부 필수] 사용자 정의 템플릿이 관련된 서열의 ID. 복합체를 모델링하는 경우에만 필요합니다. 사용자 정의 템플릿 목록을 제공하는 경우 모두 동일한 대상과 관련된 경우 단일 대상 ID를 제공할 수 있습니다. 그렇지 않으면 사용자 정의 템플릿 목록에 해당하는 대상 ID 목록을 제공해야 합니다.

#### 시각화 인수
- `--no_server`: [선택 사항] 출력 페이지에 대한 서버를 실행하지 않는 플래그 (아래 참조). 페이지는 생성됩니다. 클러스터에서 실행할 때 유용합니다.
- `--no_visuals`: [선택 사항] 출력 페이지나 PAE 플롯을 생성하지 않고 모델만 출력하는 플래그.

#### 사용자 정의 템플릿 사용법

ID가 `A`인 단백질 서열에 대해 사용자 정의 템플릿 `custom_a.pdb`를 제공하고 싶고, 템플릿에 체인 `A`와 체인 `B` 두 개의 체인이 있으며 체인 `B`를 템플릿으로 사용하고 싶다면 다음을 실행할 수 있습니다:

```bash
abcfold <input_json>  <output_dir> -abc --mmseqs2 --custom_template custom_a.pdb  --custom_template_chain B --target_id A
```

입력 서열에 여러 ID가 있고 여러 템플릿 파일이 있으며 `custom_a.pdb`에서 체인 `A`, `custom_b.pdb`에서 체인 `B`, `custom_c.pdb`에서 체인 B, 이렇게 3개의 사용자 정의 템플릿을 제공하고 싶고, 여기서 `custom_a.pdb`와 `custom_b.pdb`는 ID `A`에 해당하고 `custom_c.pdb`는 ID `B`에 해당한다면 다음을 실행할 수 있습니다:

```bash
abcfold <input_json>  <output_dir> -abc --mmseqs2 --custom_template custom_a.pdb custom_b.pdb custom_c.pdb --custom_template_chain A B B --target_id A A B
```

### 출력

ABCFold는 `<output_dir>`에 AlphaFold, Boltz 및/또는 Chai 모델을 출력하며, 결과 테이블과 유용한 [PAE 뷰어](https://gitlab.gwdg.de/general-microbiology/pae-viewer)를 포함하는 출력 페이지도 생성합니다. `--no_server` 또는 `--no_visuals` 플래그를 사용하지 않는 한 기본 브라우저에서 자동으로 열립니다.

`--no_visuals` 플래그를 사용하지 않는 한 다음을 실행하여 출력 페이지를 열 수 있습니다:

```bash
cd <output_dir>
python open_output.py
```

## 메인 페이지 예시
![main_page_example](https://raw.githubusercontent.com/rigdenlab/ABCFold/refs/heads/main/abcfold/html/static/main_page_example.png)

## PAE 뷰어 예시
![pae_viewer_example](https://raw.githubusercontent.com/rigdenlab/ABCFold/refs/heads/main/abcfold/html/static/pae_viewer_example.png)

출력 페이지는 `http://localhost:8000/index.html`에서 사용할 수 있습니다. 출력을 생성하기 위해 서버를 다시 실행해야 하는 경우 `<output_dir>`에서 `open_output.py`를 찾을 수 있습니다. 이것은 `<output_dir>`에서 실행해야 합니다.

## 추가 기능

다음은 AlphaFold3 입력 JSON 파일에 MMseqs2 MSA 및 사용자 정의 템플릿을 추가하는 스크립트입니다.

> [!WARNING]
> 이러한 스크립트는 입력 JSON 파일만 수정하며, 즉 AlphaFold3, Boltz 및 Chai-1을 실행하지 않습니다.

### MMseqs2 MSA 및 템플릿 추가

AlphaFold3 입력 JSON에 MMseqs2 MSA 및 템플릿을 추가하려면 `mmseqs2msa`를 사용할 수 있습니다:

#### 템플릿 사용

템플릿과 함께 스크립트를 실행하려면 다음 명령을 사용하십시오:

```bash
mmseqs2msa --input_json <input_json> --output_json <output_json> --templates --num_templates <num_templates>
```

- `<input_json>`: 입력 AlphaFold3 JSON 파일 경로.
- `<output_json>`: [선택 사항] 출력 JSON 파일 경로 (기본값: `<input_json_stem>`_mmseqs.json).
- `<num_templates>`: [선택 사항] 사용할 템플릿 수 (기본값: 20)
- `<mmseqs_database>`: [선택 사항] 로컬 MMSeqs2 사본에서 사용하는 데이터베이스 경로. mmseqs가 설치되어 있다면 이 플래그를 포함하면 MMseqs2를 로컬에서 실행할 수 있습니다.

> [!NOTE]
> mmseqs 데이터베이스를 설치해야 하는 경우 setup_mmseqs_databases.sh를 사용할 수 있습니다.
> 이것은 ColabFold의 MMSeqs2 데이터베이스 설정을 복제합니다.

```bash
MMSEQS_NO_INDEX=1 ./setup_mmseqs_databases.sh /path/to/db_folder
```

#### 템플릿 미사용

템플릿 없이 스크립트를 실행하려면 다음 명령을 사용하십시오:

```bash
mmseqs2msa --input_json <input_json> --output_json <output_json>
```

- `<input_json>`: 입력 AlphaFold3 JSON 파일 경로.
- `<output_json>`: [선택 사항] 출력 JSON 파일 경로 (기본값: `<input_json_stem>`_mmseqs.json).

### 사용자 정의 템플릿 추가

AlphaFold3 작업에 사용자 정의 템플릿(예: 아직 PDB에 기탁되지 않은 상동체)을 추가하고 싶을 수 있습니다. 두 가지 방법으로 수행할 수 있습니다:

#### add_custom_template.py

사용자 정의 템플릿만 추가하고 싶다면 `custom_templates`를 사용할 수 있습니다:

```bash
custom_templates --input_json <input_json> --output_json <output_json> --custom_template <custom_template> --custom_template_chain <custom_template_chain> --target_id <target_id>
```

- `<input_json>`: 입력 AlphaFold3 JSON 파일 경로.
- `<output_json>`: [선택 사항] 출력 JSON 파일 경로 (기본값: `<input_json_stem>`_custom_template.json).
- `<custom_template>`: [선택 사항] mmCIF 형식의 사용자 정의 템플릿 파일 경로 또는 사용자 정의 템플릿 목록.
- `<custom_template_chain>`: [조건부 필수] 사용자 정의 템플릿에서 사용할 체인의 체인 ID. 다중 체인 템플릿을 사용하는 경우에만 필요합니다. 사용자 정의 템플릿 목록을 제공하는 경우 사용자 정의 템플릿 체인 목록을 제공해야 합니다.
- `<target_id>`: [조건부 필수] 사용자 정의 템플릿이 관련된 서열의 ID. 복합체를 모델링하는 경우에만 필요합니다. 사용자 정의 템플릿 목록을 제공하는 경우 모두 동일한 대상과 관련된 경우 단일 대상 ID를 제공할 수 있습니다. 그렇지 않으면 사용자 정의 템플릿 목록에 해당하는 대상 ID 목록을 제공해야 합니다.

#### add_mmseqs_msa.py

사용자 정의 템플릿을 추가하고 MMseqs2 MSA/템플릿을 생성하려면 `mmseqs2msa`를 사용할 수 있습니다:

```bash
mmseqs2msa --input_json <input_json> --output_json <output_json> --templates --num_templates <num_templates> --custom_template <custom_template> --custom_template_chain <custom_template_chain> --target_id <target_id>
```

- `<input_json>`: 입력 AlphaFold3 JSON 파일 경로.
- `<output_json>`: [선택 사항] 출력 JSON 파일 경로 (기본값: `<input_json_stem>`_mmseqs.json).
- `<num_templates>`: [선택 사항] 사용할 템플릿 수 (기본값: 20)
- `<custom_template>`: [선택 사항] mmCIF 형식의 사용자 정의 템플릿 파일 경로 또는 사용자 정의 템플릿 목록.
- `<custom_template_chain>`: [조건부 필수] 사용자 정의 템플릿에서 사용할 체인의 체인 ID. 다중 체인 템플릿을 사용하는 경우에만 필요합니다. 사용자 정의 템플릿 목록을 제공하는 경우 사용자 정의 템플릿 체인 목록을 제공해야 합니다.
- `<target_id>`: [조건부 필수] 사용자 정의 템플릿이 관련된 서열의 ID. 복합체를 모델링하는 경우에만 필요합니다. 사용자 정의 템플릿 목록을 제공하는 경우 모두 동일한 대상과 관련된 경우 단일 대상 ID를 제공할 수 있습니다. 그렇지 않으면 사용자 정의 템플릿 목록에 해당하는 대상 ID 목록을 제공해야 합니다.

### 발생 가능한 문제

#### 동형 올리고머와 함께 `--target_id` 사용

아래는 이종 3량체(hetero-3-mer)의 예입니다. 동형 올리고머(homo-oligomer)를 모델링할 때 id는 목록으로 주어지며, 목록에서 식별자 중 하나를 선택해야 합니다.

```json
{
  "name": "7ZYH",
  "sequences": [
    {
      "protein": {
        "id": "A",
        "sequence": "SNAESKIKDCPWYDRGFCKHGPLCRHRHTRRVICVNYLVGFCPEGPSCKFMHPRFELPMGTTEQ"
      }
    },
    {
      "protein": {
        "id": ["B", "C"],
        "sequence": "SNAGSINGVPLLEVDLDSFEDKPWRKPGADLSDYFNYGFNEDTWKAYCEKQKRIRMGLEVIPVTSTTNK"
      }
    }
  ],
  "modelSeeds": [1],
  "dialect": "alphafold3",
  "version": 1
}
```

첫 번째 서열에 사용자 정의 템플릿을 추가하려면 `--target_id A`를 사용합니다. 두 번째 서열에 사용자 정의 템플릿을 추가하려면 `--target_id B` 또는 `--target_id C`를 사용합니다.

#### Boltz 제한 사항

Boltz에서 동일한 서열의 여러 복사본을 모델링하는 경우 입력 JSON은 다음과 같이 설정해야 합니다:

```json
{
  "name": "7ZYH",
  "sequences": [
    {
      "protein": {
        "id": ["A", "B"],
        "sequence": "SNAESKIKDCPWYDRGFCKHGPLCRHRHTRRVICVNYLVGFCPEGPSCKFMHPRFELPMGTTEQ"
      }
    }
  ],
  "modelSeeds": [1],
  "dialect": "alphafold3",
  "version": 1
}
```

동일한 서열이 아래와 같이 별개의 엔티티로 주어지면 오류가 발생합니다.

```json
{
  "name": "7ZYH",
  "sequences": [
    {
      "protein": {
        "id": "A",
        "sequence": "SNAESKIKDCPWYDRGFCKHGPLCRHRHTRRVICVNYLVGFCPEGPSCKFMHPRFELPMGTTEQ"
      }
    },
    {
      "protein": {
        "id": "B",
        "sequence": "SNAESKIKDCPWYDRGFCKHGPLCRHRHTRRVICVNYLVGFCPEGPSCKFMHPRFELPMGTTEQ"
      }
    }
  ],
  "modelSeeds": [1],
  "dialect": "alphafold3",
  "version": 1
}
```

또한 Boltz는 현재 연결된 리간드를 생성하는 기능이 없으므로 체인/리간드 간의 공유 결합이 누락됩니다.

## 기여

기여를 환영합니다! 이슈를 열거나 풀 리퀘스트를 제출해 주십시오.
