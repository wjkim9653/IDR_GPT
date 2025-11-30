import requests
import pandas as pd
from tqdm import tqdm
import io


def get_sequences_from_files(file_paths):
    # 모든 파일에서 유니크한 UniProt ID 수집
    uniprot_ids = set()

    for file_path in file_paths:
        try:
            # 텍스트 파일 읽기 (탭 구분)
            # MobiDB와 D2P2, IUPred2A 모두 3번째 컬럼이 UniProt ID임
            df = pd.read_csv(file_path, sep='\t')

            # 컬럼명이 파일마다 다를 수 있으므로 3번째 컬럼(인덱스 2)을 가져옴
            # 데이터에 헤더가 있는 경우를 가정
            ids = df.iloc[:, 2].dropna().unique()
            uniprot_ids.update(ids)
            print(f"{file_path}: {len(ids)} 개의 ID 추출됨")

        except Exception as e:
            print(f"Error reading {file_path}: {e}")

    print(f"총 {len(uniprot_ids)} 개의 고유 UniProt ID를 찾았습니다.")

    # UniProt API를 통해 FASTA 서열 다운로드
    base_url = "https://rest.uniprot.org/uniprotkb/accessions"

    # ID 리스트를 콤마로 구분된 문자열로 변환
    ids_list = list(uniprot_ids)

    # 한 번에 너무 많은 요청을 보내면 실패할 수 있으므로 500개씩 나누어 처리
    chunk_size = 500
    all_fastas = ""

    for i in tqdm(range(0, len(ids_list), chunk_size), desc="서열 다운로드 중..."):
        chunk = ids_list[i:i + chunk_size]
        params = {
            'accessions': ','.join(chunk),
            'format': 'fasta'
        }
        response = requests.get(base_url, params=params)

        if response.ok:
            all_fastas += response.text
        else:
            print(f"Error fetching chunk {i // chunk_size + 1}")

    # 결과 저장
    output_filename = "IDR_sequences.fasta"
    with open(output_filename, "w") as f:
        f.write(all_fastas)

    print(f"완료! 모든 서열이 '{output_filename}'에 저장되었습니다.")

files = [
    '/Users/wonjinkim/PycharmProjects/helpeugene/data/disordered_ids/MobiDB.txt',
    '/Users/wonjinkim/PycharmProjects/helpeugene/data/disordered_ids/DisProt.txt',
    '/Users/wonjinkim/PycharmProjects/helpeugene/data/disordered_ids/D2P2.txt',
    '/Users/wonjinkim/PycharmProjects/helpeugene/data/disordered_ids/IUPred2A.txt'
]
get_sequences_from_files(files)