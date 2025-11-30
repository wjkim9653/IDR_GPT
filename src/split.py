import pandas as pd
import re
from tqdm import tqdm

# ==========================================
# 1. 설정 (파일 경로)
# ==========================================
# 메타데이터 파일 목록
metadata_files = [
    '/Users/wonjinkim/PycharmProjects/helpeugene/data/disordered_ids/MobiDB.txt',
    '/Users/wonjinkim/PycharmProjects/helpeugene/data/disordered_ids/DisProt.txt',
    '/Users/wonjinkim/PycharmProjects/helpeugene/data/disordered_ids/D2P2.txt',
    '/Users/wonjinkim/PycharmProjects/helpeugene/data/disordered_ids/IUPred2A.txt'
]

# 아까 다운로드 받은 전체 서열 파일
fasta_file = "IDR_sequences.fasta"

# 결과 저장할 파일 이름
output_file = "final_idr_dataset.txt"


# ==========================================
# 2. FASTA 파일 파싱 함수 (전체 서열 메모리에 로드)
# ==========================================
def parse_fasta(fasta_path):
    print("FASTA 파일 로딩 중...")
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in tqdm(f, desc="parsing FASTA"):
            line = line.strip()
            if line.startswith(">"):
                # 이전 서열 저장
                if current_id:
                    sequences[current_id] = "".join(current_seq)

                # 새 헤더 파싱 (>sp|Q9Y2X7|...)
                # UniProt ID는 보통 두 번째 파이프(|) 사이에 있음
                parts = line.split('|')
                if len(parts) >= 2:
                    current_id = parts[1]  # Q9Y2X7
                else:
                    # 헤더 형식이 다를 경우 첫 단어를 ID로 사용
                    current_id = line[1:].split()[0]

                current_seq = []
            else:
                current_seq.append(line)

        # 마지막 서열 저장
        if current_id:
            sequences[current_id] = "".join(current_seq)

    print(f"-> 총 {len(sequences)} 개의 전체 단백질 서열 로드 완료.")
    return sequences


# ==========================================
# 3. 메타데이터 파싱 및 IDR 추출
# ==========================================
def extract_idrs(meta_files, full_sequences):
    idr_dataset = []
    processed_count = 0
    error_count = 0

    print("IDR 추출 작업 시작...")

    for file_path in meta_files:
        print(f"Processing {file_path}...")
        try:
            df = pd.read_csv(file_path, sep='\t')

            # 파일마다 컬럼 구조가 조금씩 다를 수 있어 처리
            # 공통적으로 3번째 컬럼(index 2)이 UniProt ID라고 가정

            # MobiDB 처리 로직
            if 'MobiDB' in file_path:
                # MobiDB 포맷: Annotation type 컬럼에 'D_WC: (1-36)' 같은 형식으로 존재
                for _, row in tqdm(df.iterrows(), desc="processing MobiDB..."):
                    uid = row.iloc[2]  # UniProt ID
                    anno = str(row.iloc[4])  # Annotation type

                    if uid not in full_sequences:
                        continue

                    # 정규표현식으로 숫자 좌표 추출 (start-end)
                    matches = re.findall(r'\((\d+)-(\d+)\)', anno)
                    for start, end in matches:
                        s, e = int(start), int(end)
                        # 파이썬 인덱스는 0부터 시작하므로 start-1
                        seq_slice = full_sequences[uid][s - 1:e]

                        # 너무 짧은 조각(예: 5개 미만)은 제외 (선택사항)
                        if len(seq_slice) >= 10:
                            idr_dataset.append(seq_slice)
                            processed_count += 1

            # D2P2 처리 로직
            elif 'D2P2' in file_path:
                # D2P2 포맷: Start(4), End(5) 컬럼이 따로 있음 (인덱스 확인 필요)
                # snippet 기준: UniProt ID(2), Start(4), End(5)
                for _, row in tqdm(df.iterrows(), desc="processing D2P2..."):
                    uid = row.iloc[2]
                    try:
                        s = int(row.iloc[4])
                        e = int(row.iloc[5])

                        if uid in full_sequences:
                            seq_slice = full_sequences[uid][s - 1:e]
                            if len(seq_slice) >= 10:
                                idr_dataset.append(seq_slice)
                                processed_count += 1
                    except:
                        pass  # 숫자가 아닌 값이 있거나 에러 발생 시 건너뜀

            # DisProt 처리 로직
            elif 'DisProt' in file_path:
                # DisProt 포맷: Sequence(6) 컬럼이 따로 있음
                for _, row in tqdm(df.iterrows(), desc="processing DisProt..."):
                    uid = row.iloc[2]
                    try:
                        # s = int(row.iloc[4])
                        # e = int(row.iloc[5])

                        if uid in full_sequences:
                            # seq_slice = full_sequences[uid][s - 1:e]
                            seq_slice = row.iloc[6]
                            if len(seq_slice) >= 10:
                                idr_dataset.append(seq_slice)
                                processed_count += 1
                    except:
                        pass  # 숫자가 아닌 값이 있거나 에러 발생 시 건너뜀

            # IUPred2A 처리 로직
            elif 'IUPred2A' in file_path:
                # 주의: IUPred 파일은 한 줄에 아미노산 하나씩(Position) 되어 있어 처리가 다름
                # 이 경우 파일 전체를 읽어 그룹화해야 하므로 메모리를 많이 씀
                # 간단하게 하기 위해: '연속된' 포지션인지 확인해야 하는데 복잡함.
                # 전략 수정: IUPred 파일은 용량이 크고 파싱이 까다로우므로
                # 이번 예시에서는 일단 제외하거나, 단순화된 로직이 필요함.
                # 사용자가 제공한 snippet만으로는 로직을 짜기 애매하여
                # 여기서는 'MobiDB'와 'D2P2' 위주로 처리하고 스킵 메시지를 남김.
                print(f"  -> {file_path}는 residue 단위 파일이라 별도 전처리가 필요하여 이번 코드에서는 건너뜁니다.")
                pass

        except Exception as e:
            print(f"  -> Error parsing {file_path}: {e}")
            error_count += 1

    return list(set(idr_dataset))  # 중복 제거 후 반환


# ==========================================
# 4. 실행 및 저장
# ==========================================
# 1) 전체 서열 로드
full_seqs = parse_fasta(fasta_file)

# 2) IDR만 잘라내기
final_idrs = extract_idrs(metadata_files, full_seqs)

# 3) 결과 저장 (ProtGPT2 학습용 포맷)
print(f"\n최종적으로 {len(final_idrs)} 개의 고유 IDR 시퀀스를 확보했습니다.")
print(f"결과를 {output_file}에 저장합니다...")

with open(output_file, "w") as f:
    for seq in final_idrs:
        # ProtGPT2 학습을 위해 각 줄 끝에 <|endoftext|> 토큰 추가 (선택사항)
        # 단순히 줄바꿈만 해도 됨. 여기서는 줄바꿈으로 처리.
        f.write(f"{seq}\n")

print("완료!")