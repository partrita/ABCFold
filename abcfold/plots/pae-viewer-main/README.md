# PAE 뷰어

이 저장소에는 [PAE 뷰어 웹서버](http://subtiwiki.uni-goettingen.de/v4/paeViewerDemo)의 주요 구성 요소에 대한 소스 코드가 포함되어 있습니다.

현재 개별 구성 요소는 안타깝게도 [SubtiWiki](http://subtiwiki.uni-goettingen.de/v4/) 프레임워크에 포함되어 있어 자체적으로 실행할 수 없습니다. HTML 템플릿도 이 프레임워크에서 제공하는 구문을 사용합니다.

MIT 라이선스로 배포됩니다.

## 오프라인 사용법
브라우저에서 PAE 뷰어를 오프라인으로 사용하고 명령줄 인터페이스로 시작하려면 [general-microbiology/pae-viewer](https://gitlab.gwdg.de/general-microbiology/pae-viewer) 저장소에서 제공하는 스크립트를 사용할 수 있습니다. 오프라인 버전을 실행하려면 Python >=3.9가 필요합니다.

### 지침
1. PAE 뷰어 프로젝트 파일 다운로드 ([pae-viewer-main.zip](https://gitlab.gwdg.de/general-microbiology/pae-viewer/-/archive/main/pae-viewer-main.zip)) 후 압축을 해제합니다.

2. 터미널을 사용하여 현재 작업 디렉터리를 프로젝트 루트(`/your/download/path/pae-viewer-main`)로 변경하고 Python으로 로컬 HTTP 서버를 시작합니다.

   ```bash
   cd /your/download/path/pae-viewer-main
   python -m http.server 8000
   ```
3. 브라우저에서 http://localhost:8000/standalone/pae-viewer.html을 열어 오프라인 버전을 시작합니다.

4. CLI를 통해 제공된 인수로 브라우저에서 PAE 뷰어를 시작하려면 프로젝트 디렉터리에서 `pae-viewer-main/standalone/pae_viewer.py` 스크립트를 실행합니다. 예시:

   ```bash
   cd /your/download/path/pae-viewer-main/standalone
   python pae_viewer.py \
     --structure pae-viewer/sample/GatA-GatB/fold_gata_gatb_model_0.cif \
     --labels 'GatA;GatB' \
     --scores pae-viewer/sample/GatA-GatB/fold_gata_gatb_full_data_0.json \
     --crosslinks pae-viewer/sample/GatA-GatB/GatA-GatB.csv
   ```

   이를 위해서는 로컬 HTTP 서버도 실행 중이어야 합니다. Python 스크립트는 모든 데이터가 포함된 HTML 세션 파일을 만들고 브라우저에서 엽니다. 기본적으로 업로드 양식을 미리 채우고 제출합니다. 세션 파일은 `pae-viewer-main/standalone` 디렉터리에 영구적으로 저장되며 언제든지 다시 방문할 수 있습니다.


## 종속성
- [Bootstrap 5.2.3](https://getbootstrap.com/)
- [Chroma 2.4.2](https://gka.github.io/chroma.js/)
- [NGL Viewer 2.0.1](http://nglviewer.org/#ngl)
- [FileSaver.js 2.0.4](https://github.com/eligrey/FileSaver.js#readme)
