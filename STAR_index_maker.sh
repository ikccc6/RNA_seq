#!/bin/bash
module load ngs/star

# kyoto university server 환경에서 적용 가능
# reference fasta 사전 다운로드 필요
# directory 설정 변경 필요

# ==========================================
# 1. 경로 설정 (슈퍼컴퓨터 환경)
# ==========================================
REF_DIR="$HOME/reference"
GENOME_DIR="$REF_DIR/genome"
STAR_INDEX_DIR="$REF_DIR/STAR_index"

echo -e "\n[ STAR Index 생성기 ]"
echo "탐색 경로: $GENOME_DIR"

# ==========================================
# 2. Genome 폴더 탐색 및 배열 저장
# ==========================================
GENOME_FOLDERS=()
if [ -d "$GENOME_DIR" ]; then
    for d in "$GENOME_DIR"/*/; do
        [ -d "$d" ] || continue
        GENOME_FOLDERS+=("$(basename "$d")")
    done
else
    echo "  ❌ [Error] $GENOME_DIR 폴더를 찾을 수 없습니다."
    exit 1
fi

# ==========================================
# 3. 폴더 번호 선택
# ==========================================

if [ ${#GENOME_FOLDERS[@]} -gt 0 ]; then
    echo "사용 가능한 Genome 목록:"
    for i in "${!GENOME_FOLDERS[@]}"; do
        num=$((i+1))
        echo "  $num) ${GENOME_FOLDERS[$i]}"
    done

    echo ""

    # 실행 인자로 번호가 전달된 경우
    if [ -n "$1" ]; then
        choice="$1"
    else
        while true; do
            read -p "▶ Index를 생성할 Genome 번호를 선택하세요: " choice

            if [[ "$choice" =~ ^[0-9]+$ ]] && \
               [ "$choice" -ge 1 ] && \
               [ "$choice" -le "${#GENOME_FOLDERS[@]}" ]; then
                break
            else
                echo "  ❌ [Error] 1에서 ${#GENOME_FOLDERS[@]} 사이의 숫자를 입력해주세요."
            fi
        done
    fi

    # 인자로 받은 값도 검증
    if [[ "$choice" =~ ^[0-9]+$ ]] && \
       [ "$choice" -ge 1 ] && \
       [ "$choice" -le "${#GENOME_FOLDERS[@]}" ]; then

        TARGET_GENOME="${GENOME_FOLDERS[$((choice-1))]}"
        echo "  ▷ 선택된 Genome: $TARGET_GENOME"

    else
        echo "  ❌ [Error] 1에서 ${#GENOME_FOLDERS[@]} 사이의 숫자를 입력해주세요."
        exit 1
    fi

else
    echo "  ❌ [Error] $GENOME_DIR 내에 하위 폴더가 없습니다."
    exit 1
fi

# ==========================================
# 4. 입력 파일 처리 및 옵션 동적 계산
# ==========================================
# FASTA 파일 스캔
TARGET_FASTA=($(ls "$GENOME_DIR/$TARGET_GENOME"/*.{fna,fasta,fa} 2>/dev/null))

if [ ${#TARGET_FASTA[@]} -eq 0 ]; then
    echo "  ❌ [Error] $TARGET_GENOME 폴더 내에 .fna 파일이 존재하지 않습니다."
    exit 1
fi

# awk를 이용한 유전체 총 염기서열 길이(bp) 계산 (헤더 '>' 라인 제외)
echo "  - 유전체 총 길이(Genome Length)를 계산 중입니다..."
GENOME_LENGTH=$(awk '!/^>/{sum+=length($0)} END{print sum}' "${TARGET_FASTA[@]}")

if [ -z "$GENOME_LENGTH" ] || [ "$GENOME_LENGTH" -eq 0 ]; then
    echo "  ❌ [Error] 유전체 길이 계산에 실패했습니다. FASTA 파일을 확인하세요."
    exit 1
fi

# --genomeSAindexNbases 계산식 적용: min(14, int(log2(Length)/2 - 1))
SA_INDEX_NBASES=$(awk -v len="$GENOME_LENGTH" 'BEGIN {
    val = (log(len)/log(2))/2 - 1;
    if (val > 14) val = 14;
    print int(val)
}')

echo "  ▷ 계산된 Genome Length: $GENOME_LENGTH bp"
echo "  ▷ 적용될 --genomeSAindexNbases: $SA_INDEX_NBASES"

# GTF/GFF 파일 스캔 (우선순위: 첫 번째 발견되는 gtf, gff, gff3 파일)
TARGET_GTF=$(ls "$GENOME_DIR/$TARGET_GENOME"/*.gtf "$GENOME_DIR/$TARGET_GENOME"/*.gff "$GENOME_DIR/$TARGET_GENOME"/*.gff3 2>/dev/null | head -n 1)

# 출력 디렉토리 생성
OUTPUT_DIR="$STAR_INDEX_DIR/$TARGET_GENOME"
mkdir -p "$OUTPUT_DIR"

read -p "▶ 사용할 스레드(Thread) 개수를 입력하세요 (예: 32, 64): " THREADS
if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
    THREADS=16
    echo "  ▷ 기본값 16 스레드로 진행합니다."
fi

# ==========================================
# 5. STAR Index 옵션 조립 및 실행
# ==========================================
echo -e "\n🚀 STAR Index 생성을 시작합니다..."
echo "  - Output DIR: $OUTPUT_DIR"

# 기본 옵션 배열 구성
STAR_ARGS=(
    STAR
    --runMode genomeGenerate
    --runThreadN "$THREADS"
    --genomeDir "$OUTPUT_DIR"
    --genomeFastaFiles "${TARGET_FASTA[@]}"
    --genomeSAindexNbases "$SA_INDEX_NBASES"
)

# GTF 파일 존재 유무에 따른 옵션 동적 추가
if [ -n "$TARGET_GTF" ]; then
    echo "  - Annotation 파일 감지: $(basename "$TARGET_GTF")"
    echo "  - ▷ [추가 옵션] --sjdbGTFfile 적용 및 --sjdbOverhang 149 (PE 150bp 기준)"
    STAR_ARGS+=(
        --sjdbGTFfile "$TARGET_GTF"
        --sjdbOverhang 149
    )
else
    echo "  - ▷ Annotation 파일 없음: 스플라이스 접합부 인덱싱을 생략합니다."
fi

# 최종 배열 명령어로 실행
"${STAR_ARGS[@]}"

echo -e "\n✅ [$TARGET_GENOME] STAR Index 생성이 완료되었습니다."
module unload ngs/star
