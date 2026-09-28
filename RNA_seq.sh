#!/bin/bash

# =====================================================================
# 0. 로깅(Logging) 함수 및 에러 방지 설정
# =====================================================================
set -e
set -u
set -o pipefail

print_start() {
    echo -e "\n\033[1;34m========================================================\033[0m"
    echo -e "🚀 [$(date +'%Y-%m-%d %H:%M:%S')] START: $1"
    echo -e "\033[1;34m========================================================\033[0m"
}

print_end() {
    echo -e "\033[1;32m✅ [$(date +'%Y-%m-%d %H:%M:%S')] DONE: $1\033[0m"
    echo -e "\033[1;34m========================================================\033[0m\n"
}

# =====================================================================
# 1. 환경 설정 및 작업 디렉토리 인식 (수정됨: 02-1, 02-2 분리)
# =====================================================================
BASE_DIR="$HOME/RNA_seq"

echo -e "\n\033[1;36m[ 분석 프로젝트 선택 ]\033[0m"
echo "${BASE_DIR} 경로의 프로젝트 폴더 목록:"

# 1. YYMMDD_name 패턴의 폴더만 스캔하여 배열에 저장
PROJECT_DIRS=()
if [ -d "$BASE_DIR" ]; then
    for d in "$BASE_DIR"/*/; do
        # 디렉토리가 없는 경우 방어 코드
        [ -d "$d" ] || continue
        
        dirname="$(basename "$d")"
        
        # 정규식: 맨 앞이 숫자 6자리로 시작하고 이어서 언더바(_)가 오는 폴더만 필터링
        if [[ "$dirname" =~ ^[0-9]{6}_ ]]; then
            PROJECT_DIRS+=("$dirname")
        fi
    done
else
    echo "  ❌ [Error] 기준 경로($BASE_DIR)를 찾을 수 없습니다."
    exit 1
fi

# 2. 번호 목록 출력 및 선택 로직
if [ ${#PROJECT_DIRS[@]} -gt 0 ]; then
    for i in "${!PROJECT_DIRS[@]}"; do
        num=$((i+1))
        echo "  $num) ${PROJECT_DIRS[$i]}"
    done
    
    echo ""
    while true; do
        read -p "▶ 작업할 프로젝트 번호를 선택하세요: " proj_choice
        
        # 입력값이 숫자인지, 범위 내에 있는지 검증
        if [[ "$proj_choice" =~ ^[0-9]+$ ]] && [ "$proj_choice" -ge 1 ] && [ "$proj_choice" -le "${#PROJECT_DIRS[@]}" ]; then
            TARGET_PROJECT="${PROJECT_DIRS[$((proj_choice-1))]}"
            echo "  ▷ 선택된 프로젝트: $TARGET_PROJECT"
            
            # 최종 프로젝트 작업 경로 변수 할당
            PROJECT_DIR="$BASE_DIR/$TARGET_PROJECT"
            break
        else
            echo "  ❌ [Error] 1에서 ${#PROJECT_DIRS[@]} 사이의 유효한 숫자를 입력해주세요."
        fi
    done
else
    echo "  ❌ [Error] $BASE_DIR 내에 YYMMDD_name 규칙을 따르는 프로젝트 폴더가 없습니다."
    exit 1
fi

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rna_seq

RAW_DIR="$PROJECT_DIR/00_raw_data"
QC_DIR="$PROJECT_DIR/01_qc"
HISAT2_DIR="$PROJECT_DIR/02-1_hisat_alignment"
STAR_DIR="$PROJECT_DIR/02-2_STAR_alignment"
ASSEMBLY_DIR="$PROJECT_DIR/03_assembly"
ANALYSIS_DIR="$PROJECT_DIR/04_analysis"
VIS_DIR="$PROJECT_DIR/05_visualization"

mkdir -p "$QC_DIR/trimmed"
mkdir -p "$QC_DIR/fastqc_trim"
mkdir -p "$QC_DIR/logs"
mkdir -p "$HISAT2_DIR"
mkdir -p "$STAR_DIR"
mkdir -p "$ASSEMBLY_DIR"
mkdir -p "$ANALYSIS_DIR"
mkdir -p "$VIS_DIR"

print_end "디렉토리 준비 완료"

# =====================================================================
# 2. 샘플 리스트 자동 추출
# =====================================================================
print_start "샘플 목록 추출"

SAMPLE_LIST=()

for FILE in "$RAW_DIR"/*.f*q*; do
    [ -e "$FILE" ] || continue 
    
    FILENAME=$(basename "$FILE")
    SAMPLE_ID=$(echo "$FILENAME" | sed -E 's/(_L[0-9]+)?([._][R]?[12])?\.f(ast)?q(\.gz)?//')
    
    SAMPLE_LIST+=("$SAMPLE_ID")
done

UNIQUE_SAMPLES=($(printf "%s\n" "${SAMPLE_LIST[@]}" | sort -u))

echo "총 ${#UNIQUE_SAMPLES[@]}개의 샘플이 확인되었습니다:"
printf " - %s\n" "${UNIQUE_SAMPLES[@]}"

print_end "샘플 목록 추출 완료"

# =====================================================================
# 3. QC (Quality Control) - fastp
# =====================================================================
print_start "Quality Control [fastp]"

for SAMPLE in "${UNIQUE_SAMPLES[@]}"; do
    echo "  ▶ Processing sample : ${SAMPLE}"

    LEGACY_QC_FILE="$QC_DIR/${SAMPLE}_filtered.fastq"
    if [ -s "$LEGACY_QC_FILE" ]; then
        echo "  ⏭️ [Skip] $SAMPLE : 기존 QC 완료 파일이 확인되었습니다. (${LEGACY_QC_FILE})"
        ln -sf "$LEGACY_QC_FILE" "$QC_DIR/trimmed/${SAMPLE}_trimmed.fastq"
        continue
    fi

    EXPECTED_PE_R1="$QC_DIR/trimmed/${SAMPLE}_1_P.fastq.gz"
    EXPECTED_SE="$QC_DIR/trimmed/${SAMPLE}_trimmed.fastq.gz"

    if [ -s "$EXPECTED_PE_R1" ] || [ -s "$EXPECTED_SE" ]; then
        echo "  ⏭️ [Skip] $SAMPLE : 이미 Trimming이 완료되었습니다."
        continue
    fi

    R2_FILE=$(ls "$RAW_DIR"/${SAMPLE}* 2>/dev/null | grep -E '(_2\.|_R2\.|_2_|_R2_|\.2\.f)' | head -n 1)
    FASTP_LOG="$QC_DIR/logs/${SAMPLE}_fastp.log"

    if [ -z "$R2_FILE" ]; then
        echo "  ▷ [SE] Single-end 데이터로 식별되었습니다."
        R1_FILE=$(ls "$RAW_DIR"/${SAMPLE}* 2>/dev/null | head -n 1)
        
        fastp \
            -i "${R1_FILE}" \
            -o "${EXPECTED_SE}" \
            -h "$QC_DIR/fastqc_trim/${SAMPLE}_report.html" \
            -j "$QC_DIR/fastqc_trim/${SAMPLE}_report.json" \
            -w 16 2> "$FASTP_LOG"
    else
        echo "  ▷ [PE] Paired-end 데이터로 식별되었습니다."
        R1_FILE=$(ls "$RAW_DIR"/${SAMPLE}* 2>/dev/null | grep -E '(_1\.|_R1\.|_1_|_R1_|\.1\.f)' | head -n 1)
        EXPECTED_PE_R2="$QC_DIR/trimmed/${SAMPLE}_2_P.fastq.gz"

        fastp --detect_adapter_for_pe \
            -i "${R1_FILE}" \
            -I "${R2_FILE}" \
            -o "${EXPECTED_PE_R1}" \
            -O "${EXPECTED_PE_R2}" \
            -h "$QC_DIR/fastqc_trim/${SAMPLE}_report.html" \
            -j "$QC_DIR/fastqc_trim/${SAMPLE}_report.json" \
            -w 16 2> "$FASTP_LOG"
    fi

    echo "  ▷ [${SAMPLE}] Trimming completed"
done

print_end "Quality Control [fastp]"

# =====================================================================
# 4. 분석 단계 및 파이프라인 옵션 선택 (수정됨: HISAT2/STAR 분리)
# =====================================================================
echo -e "\n[ 분석 파이프라인 설정 ]"

read -p "▶ Alignment (HISAT2) 분석을 진행하시겠습니까? (y/n): " RUN_HISAT2
read -p "▶ Alignment (STAR) 분석을 진행하시겠습니까? (y/n): " RUN_STAR
read -p "▶ Assembly (SPAdes) 분석을 진행하시겠습니까? (y/n): " RUN_ASSEMBLY

SPADES_BASIC_OPT=""
SPADES_ADD_OPT=""
AVAIL_THREADS=$(nproc)

if [[ "$RUN_ASSEMBLY" =~ ^[Yy]$ ]]; then
    echo -e "\n[ SPAdes Basic Option 선택 ]"
    echo "  1) --rna       (RNA-Seq 데이터 일반)"
    echo "  2) --rnaviral  (RNA-Seq 데이터 기반 바이러스 탐지)"
    echo "  3) --metaviral (메타게놈 데이터 기반 바이러스 탐지)"
    echo "  4) --isolate   (High-coverage isolate 데이터)"
    echo "  5) --meta      (Metagenomic 데이터)"
    echo "  6) 사용 안 함  (직접 입력)"
    read -p "▶ 번호를 선택하세요 (1-6): " SPADES_MODE_NUM

    case $SPADES_MODE_NUM in
        1) SPADES_BASIC_OPT="--rna" ;;
        2) SPADES_BASIC_OPT="--rnaviral" ;;
        3) SPADES_BASIC_OPT="--metaviral" ;;
        4) SPADES_BASIC_OPT="--isolate" ;;
        5) SPADES_BASIC_OPT="--meta" ;;
        *) SPADES_BASIC_OPT="" ;;
    esac

    echo -e "\n[ SPAdes Pipeline/Advanced Option 추가 ]"
    echo "  * 시스템 스레드 수(-t $AVAIL_THREADS)가 자동으로 적용됩니다."
    echo "  * 예시: --only-assembler -k 21,33,55"
    
    while true; do
        read -p "▶ 추가 옵션을 입력하세요 (없으면 엔터): " SPADES_ADD_OPT

        if [[ "$SPADES_BASIC_OPT" == "--rna" || "$SPADES_BASIC_OPT" == "--rnaviral" || "$SPADES_BASIC_OPT" == "--isolate" ]]; then
            if [[ "$SPADES_ADD_OPT" == *"--only-error-correction"* || "$SPADES_ADD_OPT" == *"--careful"* ]]; then
                echo -e "  ⚠️ [경고] 선택하신 옵션($SPADES_BASIC_OPT)은 --only-error-correction 또는 --careful과 함께 사용할 수 없습니다."
                echo -e "  옵션을 다시 입력해 주세요.\n"
                continue
            fi
        fi
        break
    done
fi

# =====================================================================
# 5. Assembly (SPAdes)
# =====================================================================
if [[ "$RUN_ASSEMBLY" =~ ^[Yy]$ ]]; then
    print_start "Assembly [SPAdes]"

    mkdir -p "$ASSEMBLY_DIR/logs"

    for SAMPLE in "${UNIQUE_SAMPLES[@]}"; do
        echo "  ▶ Processing sample : ${SAMPLE}"
        
        SAMPLE_ASSEMBLY_DIR="$ASSEMBLY_DIR/$SAMPLE"
        EXPECTED_CONTIGS="$SAMPLE_ASSEMBLY_DIR/contigs.fasta"

        if [ -s "$EXPECTED_CONTIGS" ]; then
            echo "  ⏭️ [Skip] $SAMPLE : 이미 Assembly(contigs.fasta)가 완료되었습니다."
            continue
        fi

        TRIM_PE_R1="$QC_DIR/trimmed/${SAMPLE}_1_P.fastq.gz"
        TRIM_PE_R2="$QC_DIR/trimmed/${SAMPLE}_2_P.fastq.gz"
        TRIM_SE="$QC_DIR/trimmed/${SAMPLE}_trimmed.fastq.gz"

        mkdir -p "$SAMPLE_ASSEMBLY_DIR"
        SPADES_LOG="$ASSEMBLY_DIR/logs/${SAMPLE}_spades.log"

        if [ -s "$TRIM_PE_R1" ] && [ -s "$TRIM_PE_R2" ]; then
            spades.py $SPADES_BASIC_OPT \
                -1 "$TRIM_PE_R1" \
                -2 "$TRIM_PE_R2" \
                -o "$SAMPLE_ASSEMBLY_DIR" \
                -t "$AVAIL_THREADS" $SPADES_ADD_OPT > "$SPADES_LOG" 2>&1

        elif [ -s "$TRIM_SE" ]; then
            spades.py $SPADES_BASIC_OPT \
                -s "$TRIM_SE" \
                -o "$SAMPLE_ASSEMBLY_DIR" \
                -t "$AVAIL_THREADS" $SPADES_ADD_OPT > "$SPADES_LOG" 2>&1
        else
            echo "  ❌ [Error] $SAMPLE 의 QC 완료 파일을 찾을 수 없어 Assembly를 건너뜁니다."
        fi
        
        echo "  ▷ [${SAMPLE}] Assembly completed"
    done

    print_end "Assembly [SPAdes]"
else
    echo -e "\n⏭️ Assembly 단계를 건너뜁니다 (사용자 선택)."
fi

# =====================================================================
# 6. Alignment (HISAT2)
# =====================================================================
if [[ "$RUN_HISAT2" =~ ^[Yy]$ ]]; then
    print_start "Alignment [HISAT2]"

    mkdir -p "$HISAT2_DIR/logs"
    INDEX_BASE_DIR="$HOME/RNA_seq/reference/hisat_index"

    echo -e "\n[ HISAT2 Index 경로 설정 ]"
    
    # 1. 인덱스 디렉토리 스캔 및 배열 저장
    AVAILABLE_INDICES=()
    if [ -d "$INDEX_BASE_DIR" ]; then
        for d in "$INDEX_BASE_DIR"/*/; do
            if [ -d "$d" ]; then
                AVAILABLE_INDICES+=("$(basename "$d")")
            fi
        done
    else
        echo "  ❌ [Error] 인덱스 기준 폴더($INDEX_BASE_DIR)를 찾을 수 없습니다."
    fi

    # 2. 인덱스 목록 출력 및 매핑 대상 번호 선택
    if [ ${#AVAILABLE_INDICES[@]} -gt 0 ]; then
        echo "사용 가능한 HISAT2 Index 목록:"
        for i in "${!AVAILABLE_INDICES[@]}"; do
            num=$((i+1))
            echo "  $num) ${AVAILABLE_INDICES[$i]}"
        done
        
        echo ""
        while true; do
            read -p "▶ 매핑할 대상(Index) 번호를 선택하세요: " choice
            if [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -le "${#AVAILABLE_INDICES[@]}" ]; then
                INDEX_NAME="${AVAILABLE_INDICES[$((choice-1))]}"
                HISAT2_INDEX="$INDEX_BASE_DIR/$INDEX_NAME/$INDEX_NAME"
                echo "  ▷ 선택된 Index: $INDEX_NAME"
                break
            else
                echo "  ❌ [Error] 1에서 ${#AVAILABLE_INDICES[@]} 사이의 유효한 숫자를 입력해주세요."
            fi
        done

        if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
            echo "  ❌ [Error] 인덱스 파일(${HISAT2_INDEX}.1.ht2)이 존재하지 않습니다. 단계를 건너뜁니다."
        else
            echo -e "\n[ 입력 데이터 선택 ]"
            echo "  1) Trimmed FASTQ (기본 1차 정렬)"
            echo "  2) 이전 Alignment에서 추출된 Unmapped BAM (2차 정렬용)"
            read -p "▶ 번호를 선택하세요 (1 또는 2): " INPUT_TYPE
            
            PREV_INDEX=""
            if [ "$INPUT_TYPE" == "2" ]; then
                # 3. 2차 정렬용 이전 인덱스 번호 선택 로직 적용
                while true; do
                    read -p "▶ 이전에 매핑 실패했던 Index 번호를 선택하세요 (위 목록 참조): " prev_choice
                    if [[ "$prev_choice" =~ ^[0-9]+$ ]] && [ "$prev_choice" -ge 1 ] && [ "$prev_choice" -le "${#AVAILABLE_INDICES[@]}" ]; then
                        PREV_INDEX="${AVAILABLE_INDICES[$((prev_choice-1))]}"
                        echo "  ▷ 선택된 이전 Index: $PREV_INDEX"
                        break
                    else
                        echo "  ❌ [Error] 1에서 ${#AVAILABLE_INDICES[@]} 사이의 유효한 숫자를 입력해주세요."
                    fi
                done
            fi

            read -p "▶ 이번 Alignment 결과에서도 Unmapped Read를 별도로 추출하시겠습니까? (y/n): " EXTRACT_UNMAPPED

            SAMTOOLS_THREADS=4
            HISAT2_THREADS=$(( AVAIL_THREADS > 4 ? AVAIL_THREADS - 4 : 1 ))

            # 4. 샘플 반복 처리 구문
            for SAMPLE in "${UNIQUE_SAMPLES[@]}"; do
                echo "  ▶ Processing sample : ${SAMPLE}"
                
                EXPECTED_BAM="$HISAT2_DIR/${SAMPLE}_${INDEX_NAME}_aligned.bam"
                UNMAPPED_BAM="$HISAT2_DIR/${SAMPLE}_${INDEX_NAME}_unmapped.bam"
                
                ALIGN_LOG="$HISAT2_DIR/logs/${SAMPLE}_${INDEX_NAME}_hisat2.log"
                SUMMARY_FILE="$HISAT2_DIR/logs/${SAMPLE}_${INDEX_NAME}_summary.txt"

                if [ -s "$EXPECTED_BAM" ]; then
                    echo "  ⏭️ [Skip] $SAMPLE : 이미 해당 인덱스(${INDEX_NAME})에 대한 Alignment가 완료되었습니다."
                    continue
                fi

                IS_PE=0
                TRIM_PE_R1="$QC_DIR/trimmed/${SAMPLE}_1_P.fastq.gz"
                TRIM_PE_R2="$QC_DIR/trimmed/${SAMPLE}_2_P.fastq.gz"
                TRIM_SE="$QC_DIR/trimmed/${SAMPLE}_trimmed.fastq.gz"
                LEGACY_SE="$QC_DIR/trimmed/${SAMPLE}_trimmed.fastq"

                if [ -s "$TRIM_PE_R1" ] && [ -s "$TRIM_PE_R2" ]; then IS_PE=1; fi

                if [ "$INPUT_TYPE" == "1" ]; then
                    if [ "$IS_PE" -eq 1 ]; then
                        hisat2 -p "$HISAT2_THREADS" --dta --summary-file "$SUMMARY_FILE" --new-summary \
                               -x "$HISAT2_INDEX" -1 "$TRIM_PE_R1" -2 "$TRIM_PE_R2" 2>> "$ALIGN_LOG" | \
                        samtools sort -@ "$SAMTOOLS_THREADS" -o "$EXPECTED_BAM"
                    elif [ -s "$TRIM_SE" ]; then
                        hisat2 -p "$HISAT2_THREADS" --dta --summary-file "$SUMMARY_FILE" --new-summary \
                               -x "$HISAT2_INDEX" -U "$TRIM_SE" 2>> "$ALIGN_LOG" | \
                        samtools sort -@ "$SAMTOOLS_THREADS" -o "$EXPECTED_BAM"
                    elif [ -s "$LEGACY_SE" ]; then
                        hisat2 -p "$HISAT2_THREADS" --dta --summary-file "$SUMMARY_FILE" --new-summary \
                               -x "$HISAT2_INDEX" -U "$LEGACY_SE" 2>> "$ALIGN_LOG" | \
                        samtools sort -@ "$SAMTOOLS_THREADS" -o "$EXPECTED_BAM"
                    else
                        echo "  ❌ [Error] $SAMPLE 의 QC 완료 파일을 찾을 수 없습니다."
                        continue
                    fi

                elif [ "$INPUT_TYPE" == "2" ]; then
                    PREV_BAM="$HISAT2_DIR/${SAMPLE}_${PREV_INDEX}_unmapped.bam"
                    if [ ! -s "$PREV_BAM" ]; then
                        echo "  ❌ [Error] $PREV_BAM 파일을 찾을 수 없습니다. 건너뜁니다."
                        continue
                    fi

                    TMP_R1="$HISAT2_DIR/${SAMPLE}_tmp_R1.fastq"
                    TMP_R2="$HISAT2_DIR/${SAMPLE}_tmp_R2.fastq"
                    TMP_SE="$HISAT2_DIR/${SAMPLE}_tmp_SE.fastq"

                    echo "    - BAM to FASTQ 변환 중..."
                    if [ "$IS_PE" -eq 1 ]; then
                        # Pair 동기화(collate) 및 이름 무결성(-n) 유지
                        samtools collate -@ "$SAMTOOLS_THREADS" -u -O "$PREV_BAM" | \
                        samtools fastq -@ "$SAMTOOLS_THREADS" -n -1 "$TMP_R1" -2 "$TMP_R2" -0 /dev/null -s /dev/null -
                        
                        hisat2 -p "$HISAT2_THREADS" --dta --summary-file "$SUMMARY_FILE" --new-summary \
                               -x "$HISAT2_INDEX" -1 "$TMP_R1" -2 "$TMP_R2" 2>> "$ALIGN_LOG" | \
                        samtools sort -@ "$SAMTOOLS_THREADS" -o "$EXPECTED_BAM"
                        
                        rm -f "$TMP_R1" "$TMP_R2"
                    else
                        # Single-end 리드 이름 무결성(-n) 유지
                        samtools fastq -@ "$SAMTOOLS_THREADS" -n -0 "$TMP_SE" "$PREV_BAM"
                        
                        hisat2 -p "$HISAT2_THREADS" --dta --summary-file "$SUMMARY_FILE" --new-summary \
                               -x "$HISAT2_INDEX" -U "$TMP_SE" 2>> "$ALIGN_LOG" | \
                        samtools sort -@ "$SAMTOOLS_THREADS" -o "$EXPECTED_BAM"
                        
                        rm -f "$TMP_SE"
                    fi
                fi
                
                samtools index -@ "$SAMTOOLS_THREADS" "$EXPECTED_BAM"
                echo "  ▷ [${SAMPLE}] Alignment to ${INDEX_NAME} completed"
                
                if [[ "$EXTRACT_UNMAPPED" =~ ^[Yy]$ ]]; then
                    echo "    - Extracting unmapped reads..."
                    if [ "$IS_PE" -eq 1 ]; then
                        UNMAP_FLAG=12
                    else
                        UNMAP_FLAG=4
                    fi

                    samtools view -@ "$SAMTOOLS_THREADS" -b -f "$UNMAP_FLAG" -o "$UNMAPPED_BAM" "$EXPECTED_BAM"
                    samtools index -@ "$SAMTOOLS_THREADS" "$UNMAPPED_BAM"
                    echo "  ▷ [${SAMPLE}] Unmapped BAM created: ${SAMPLE}_${INDEX_NAME}_unmapped.bam"
                fi
            done
        fi
    else
        echo "  ❌ [Error] 사용 가능한 인덱스가 없으므로 단계를 건너뜁니다."
    fi

    print_end "Alignment [HISAT2]"
else
    echo -e "\n⏭️ Alignment (HISAT2) 단계를 건너뜁니다 (사용자 선택)."
fi

# =====================================================================
# 7. Alignment (STAR) - 능동적 1차/2차 정렬 흐름 적용
# =====================================================================
if [[ "$RUN_STAR" =~ ^[Yy]$ ]]; then
    print_start "Alignment [STAR]"

    mkdir -p "$STAR_DIR/logs"
    STAR_INDEX_BASE_DIR="$HOME/RNA_seq/reference/STAR_index"

    echo -e "\n[ STAR Index 경로 설정 ]"
    
    # 1. 인덱스 디렉토리 스캔 및 배열 저장
    AVAILABLE_STAR_INDICES=()
    if [ -d "$STAR_INDEX_BASE_DIR" ]; then
        for d in "$STAR_INDEX_BASE_DIR"/*/; do
            [ -d "$d" ] || continue
            AVAILABLE_STAR_INDICES+=("$(basename "$d")")
        done
    else
        echo "  ❌ [Error] 인덱스 기준 폴더($STAR_INDEX_BASE_DIR)를 찾을 수 없습니다."
    fi

    # 2. 인덱스 목록 출력 및 매핑 대상 번호 선택
    if [ ${#AVAILABLE_STAR_INDICES[@]} -gt 0 ]; then
        echo "사용 가능한 STAR Index 목록:"
        for i in "${!AVAILABLE_STAR_INDICES[@]}"; do
            num=$((i+1))
            echo "  $num) ${AVAILABLE_STAR_INDICES[$i]}"
        done
        
        echo ""
        while true; do
            read -p "▶ 매핑할 대상(Index) 번호를 선택하세요: " choice
            if [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -le "${#AVAILABLE_STAR_INDICES[@]}" ]; then
                INDEX_NAME="${AVAILABLE_STAR_INDICES[$((choice-1))]}"
                STAR_INDEX="$STAR_INDEX_BASE_DIR/$INDEX_NAME"
                echo "  ▷ 선택된 Index: $INDEX_NAME"
                break
            else
                echo "  ❌ [Error] 1에서 ${#AVAILABLE_STAR_INDICES[@]} 사이의 유효한 숫자를 입력해주세요."
            fi
        done

        echo -e "\n[ 입력 데이터 선택 ]"
        echo "  1) Trimmed FASTQ (기본 1차 정렬)"
        echo "  2) 이전 Alignment에서 추출된 Unmapped BAM (2차 정렬용)"
        read -p "▶ 번호를 선택하세요 (1 또는 2): " INPUT_TYPE
        
        PREV_INDEX=""
        if [ "$INPUT_TYPE" == "2" ]; then
            # 3. 2차 정렬용 이전 인덱스 번호 선택 로직
            while true; do
                read -p "▶ 이전에 매핑 실패했던 Index 번호를 선택하세요 (위 목록 참조): " prev_choice
                if [[ "$prev_choice" =~ ^[0-9]+$ ]] && [ "$prev_choice" -ge 1 ] && [ "$prev_choice" -le "${#AVAILABLE_STAR_INDICES[@]}" ]; then
                    PREV_INDEX="${AVAILABLE_STAR_INDICES[$((prev_choice-1))]}"
                    echo "  ▷ 선택된 이전 Index: $PREV_INDEX"
                    break
                else
                    echo "  ❌ [Error] 1에서 ${#AVAILABLE_STAR_INDICES[@]} 사이의 유효한 숫자를 입력해주세요."
                fi
            done
        fi

        read -p "▶ 이번 Alignment 결과에서도 Unmapped Read를 별도로 추출하시겠습니까? (y/n): " EXTRACT_UNMAPPED

        SAMTOOLS_THREADS=4
        STAR_THREADS=$(( AVAIL_THREADS > 4 ? AVAIL_THREADS - 4 : 1 ))
        BAM_SORT_RAM=31474836480

        # 4. 샘플 반복 처리 구문
        for SAMPLE in "${UNIQUE_SAMPLES[@]}"; do
            echo -e "\n  ▶ Processing sample : ${SAMPLE}"
            
            EXPECTED_BAM="$STAR_DIR/${SAMPLE}_${INDEX_NAME}_aligned.bam"
            UNMAPPED_BAM="$STAR_DIR/${SAMPLE}_${INDEX_NAME}_unmapped.bam"
            
            STAR_PREFIX="$STAR_DIR/${SAMPLE}_${INDEX_NAME}_"
            STAR_LOG="$STAR_DIR/logs/${SAMPLE}_${INDEX_NAME}_STAR_Log.final.out"

            if [ -s "$EXPECTED_BAM" ]; then
                echo "  ⏭️ [Skip] $SAMPLE : 이미 해당 인덱스(${INDEX_NAME})에 대한 Alignment가 완료되었습니다."
                continue
            fi

            IS_PE=0
            TRIM_PE_R1="$QC_DIR/trimmed/${SAMPLE}_1_P.fastq.gz"
            TRIM_PE_R2="$QC_DIR/trimmed/${SAMPLE}_2_P.fastq.gz"
            TRIM_SE="$QC_DIR/trimmed/${SAMPLE}_trimmed.fastq.gz"
            LEGACY_SE="$QC_DIR/trimmed/${SAMPLE}_filtered.fastq"

            if [ -s "$TRIM_PE_R1" ] && [ -s "$TRIM_PE_R2" ]; then IS_PE=1; fi

            # 공통 STAR 기본 명령어 구성
            STAR_CMD=(
                STAR
                --runThreadN "$STAR_THREADS"
                --genomeDir "$STAR_INDEX"
                --outSAMtype BAM SortedByCoordinate
                --limitBAMsortRAM "$BAM_SORT_RAM"
                --outFileNamePrefix "$STAR_PREFIX"
                --quantMode GeneCounts
            )

            if [ "$INPUT_TYPE" == "1" ]; then
                if [ "$IS_PE" -eq 1 ]; then
                    STAR_CMD+=(--readFilesIn "$TRIM_PE_R1" "$TRIM_PE_R2" --readFilesCommand zcat)
                elif [ -s "$TRIM_SE" ]; then
                    STAR_CMD+=(--readFilesIn "$TRIM_SE" --readFilesCommand zcat)
                elif [ -s "$LEGACY_SE" ]; then
                    STAR_CMD+=(--readFilesIn "$LEGACY_SE")
                else
                    echo "  ❌ [Error] $SAMPLE 의 QC 완료 파일을 찾을 수 없습니다."
                    continue
                fi

                echo "    - 1차 매핑 진행 중 ($INDEX_NAME)..."
                "${STAR_CMD[@]}" > /dev/null 2>&1

            elif [ "$INPUT_TYPE" == "2" ]; then
                PREV_BAM="$STAR_DIR/${SAMPLE}_${PREV_INDEX}_unmapped.bam"
                if [ ! -s "$PREV_BAM" ]; then
                    echo "  ❌ [Error] $PREV_BAM 파일을 찾을 수 없습니다. 건너뜁니다."
                    continue
                fi

                TMP_R1="$STAR_DIR/${SAMPLE}_tmp_R1.fastq"
                TMP_R2="$STAR_DIR/${SAMPLE}_tmp_R2.fastq"
                TMP_SE="$STAR_DIR/${SAMPLE}_tmp_SE.fastq"

                echo "    - BAM to FASTQ 변환 중..."
                if [ "$IS_PE" -eq 1 ]; then
                    # Pair 동기화(collate) 및 이름 무결성(-n) 유지
                    samtools collate -@ "$SAMTOOLS_THREADS" -u -O "$PREV_BAM" | \
                    samtools fastq -@ "$SAMTOOLS_THREADS" -n -1 "$TMP_R1" -2 "$TMP_R2" -0 /dev/null -s /dev/null -
                    
                    STAR_CMD+=(--readFilesIn "$TMP_R1" "$TMP_R2")
                    echo "    - 2차 매핑 진행 중 ($INDEX_NAME)..."
                    "${STAR_CMD[@]}" > /dev/null 2>&1
                    
                    rm -f "$TMP_R1" "$TMP_R2"
                else
                    # Single-end 리드 이름 무결성(-n) 유지
                    samtools fastq -@ "$SAMTOOLS_THREADS" -n -0 "$TMP_SE" "$PREV_BAM"
                    
                    STAR_CMD+=(--readFilesIn "$TMP_SE")
                    echo "    - 2차 매핑 진행 중 ($INDEX_NAME)..."
                    "${STAR_CMD[@]}" > /dev/null 2>&1
                    
                    rm -f "$TMP_SE"
                fi
            fi
            
            # 생성된 파일 정리 및 이름 변경
            STAR_OUTPUT_BAM="${STAR_PREFIX}Aligned.sortedByCoord.out.bam"
            if [ -f "$STAR_OUTPUT_BAM" ]; then
                mv "$STAR_OUTPUT_BAM" "$EXPECTED_BAM"
                mv "${STAR_PREFIX}Log.final.out" "$STAR_LOG"
                rm -f "${STAR_PREFIX}Log.out" "${STAR_PREFIX}Log.progress.out" "${STAR_PREFIX}SJ.out.tab"
                
                samtools index -@ "$SAMTOOLS_THREADS" "$EXPECTED_BAM"
                echo "  ▷ [${SAMPLE}] Alignment to ${INDEX_NAME} completed"
            else
                echo "  ❌ [Error] 매핑 실패."
                continue
            fi
            
            # Unmapped 리드 추출 (BAM 형태)
            if [[ "$EXTRACT_UNMAPPED" =~ ^[Yy]$ ]]; then
                echo "    - Extracting unmapped reads..."
                if [ "$IS_PE" -eq 1 ]; then
                    UNMAP_FLAG=12
                else
                    UNMAP_FLAG=4
                fi

                samtools view -@ "$SAMTOOLS_THREADS" -b -f "$UNMAP_FLAG" -o "$UNMAPPED_BAM" "$EXPECTED_BAM"
                samtools index -@ "$SAMTOOLS_THREADS" "$UNMAPPED_BAM"
                echo "  ▷ [${SAMPLE}] Unmapped BAM created: ${SAMPLE}_${INDEX_NAME}_unmapped.bam"
            fi
        done
    else
        echo "  ❌ [Error] 사용 가능한 인덱스가 없으므로 단계를 건너뜁니다."
    fi

    print_end "Alignment [STAR]"
else
    echo -e "\n⏭️ Alignment (STAR) 단계를 건너뜁니다 (사용자 선택)."
fi
