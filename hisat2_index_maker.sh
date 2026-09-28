#!/bin/bash

# =====================================================================
# HISAT2 레퍼런스 인덱스 설정 및 생성 (index_maker.sh)
# 여러 FASTA 파일을 쉼표(,)로 연결하여 HISAT2에 직접 입력
# =====================================================================
REF_GENOME_BASE="$HOME/RNA_seq/reference/genome"
REF_INDEX_BASE="$HOME/RNA_seq/reference/hisat_index"

echo -e "\n[ HISAT2 Reference Index 설정 ]"
echo "사용 가능한 레퍼런스 게놈 목록:"
ls -1 "$REF_GENOME_BASE"

read -p "▶ 사용할 레퍼런스 폴더명을 입력하세요 (예: WSN 또는 canfam3.1): " GENOME_NAME

SELECTED_GENOME_DIR="$REF_GENOME_BASE/$GENOME_NAME"
INDEX_OUT_DIR="$REF_INDEX_BASE/$GENOME_NAME"
INDEX_PREFIX="$INDEX_OUT_DIR/$GENOME_NAME"

if [ ! -d "$SELECTED_GENOME_DIR" ]; then
    echo "❌ [Error] $SELECTED_GENOME_DIR 디렉토리를 찾을 수 없습니다."
    exit 1
fi

mkdir -p "$INDEX_OUT_DIR"

if [ -s "${INDEX_PREFIX}.1.ht2" ]; then
    echo "⏭️ [Skip] $GENOME_NAME 의 HISAT2 인덱스가 이미 존재합니다. ($INDEX_PREFIX)"
else
    echo "▶ $GENOME_NAME 의 HISAT2 인덱스를 생성합니다..."
    
    shopt -s nullglob
    FASTA_FILES=(
        "$SELECTED_GENOME_DIR"/*.fasta "$SELECTED_GENOME_DIR"/*.fasta.gz
        "$SELECTED_GENOME_DIR"/*.fa "$SELECTED_GENOME_DIR"/*.fa.gz
        "$SELECTED_GENOME_DIR"/*.fna "$SELECTED_GENOME_DIR"/*.fna.gz
    )
    shopt -u nullglob

    if [ ${#FASTA_FILES[@]} -eq 0 ]; then
        echo "❌ [Error] $SELECTED_GENOME_DIR 내에 유효한 레퍼런스 파일이 존재하지 않습니다."
        exit 1
    fi
    
    # 1. 입력 파일 리스트 구성 (.gz 파일은 임시 압축 해제)
    PROCESS_FILES=()
    TEMP_DECOMPRESS_DIR="$INDEX_OUT_DIR/temp_unzip"
    mkdir -p "$TEMP_DECOMPRESS_DIR"
    
    for f in "${FASTA_FILES[@]}"; do
        if [[ "$f" == *.gz ]]; then
            BASENAME=$(basename "$f" .gz)
            gzip -dc "$f" > "$TEMP_DECOMPRESS_DIR/$BASENAME"
            PROCESS_FILES+=("$TEMP_DECOMPRESS_DIR/$BASENAME")
        else
            PROCESS_FILES+=("$f")
        fi
    done

    # 2. 파일 배열을 쉼표(,)로 연결된 하나의 문자열로 변환
    COMMA_SEPARATED_LIST=$(IFS=, ; echo "${PROCESS_FILES[*]}")
    
    echo "  - 입력 파일 리스트 연결 완료 (총 ${#PROCESS_FILES[@]}개 파일)"
    
    # 3. HISAT2 인덱스 빌드 실행 (쉼표로 연결된 리스트 직접 입력)
    hisat2-build -p $(nproc) "$COMMA_SEPARATED_LIST" "$INDEX_PREFIX"
    BUILD_STATUS=$?
    
    # 4. 임시 압축 해제 폴더 정리
    rm -rf "$TEMP_DECOMPRESS_DIR"
    
    # 5. 생성 결과 확인
    if [ $BUILD_STATUS -eq 0 ] && [ -s "${INDEX_PREFIX}.1.ht2" ]; then
        echo "▷ $GENOME_NAME 인덱스 생성 완료"
    else
        echo "❌ [Error] 인덱스 생성 중 오류가 발생했거나 파일이 정상적으로 만들어지지 않았습니다."
        exit 1
    fi
fi

