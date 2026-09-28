#!/bin/bash
set -e
set -u
set -o pipefail

# =====================================================================
# 1. 환경 설정 및 프로젝트 선택
# =====================================================================
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rna_seq

BASE_DIR="$HOME/RNA_seq"
REF_DIR="$BASE_DIR/reference/genome"

echo -e "\n[ 분석 프로젝트 선택 ]"
PROJECT_DIRS=()
if [ -d "$BASE_DIR" ]; then
    for d in "$BASE_DIR"/*/; do
        [ -d "$d" ] || continue
        dirname="$(basename "$d")"
        if [[ "$dirname" =~ ^[0-9]{6}_ ]]; then
            PROJECT_DIRS+=("$dirname")
        fi
    done
else
    echo "❌ [Error] 기준 경로($BASE_DIR)를 찾을 수 없습니다."
    exit 1
fi

if [ ${#PROJECT_DIRS[@]} -gt 0 ]; then
    for i in "${!PROJECT_DIRS[@]}"; do
        num=$((i+1))
        echo "  $num) ${PROJECT_DIRS[$i]}"
    done
    echo ""
    while true; do
        read -p "▶ 작업할 프로젝트 번호를 선택하세요: " proj_choice
        if [[ "$proj_choice" =~ ^[0-9]+$ ]] && [ "$proj_choice" -ge 1 ] && [ "$proj_choice" -le "${#PROJECT_DIRS[@]}" ]; then
            TARGET_PROJECT="${PROJECT_DIRS[$((proj_choice-1))]}"
            PROJECT_DIR="$BASE_DIR/$TARGET_PROJECT"
            echo "  ▷ 선택된 프로젝트: $TARGET_PROJECT"
            break
        else
            echo "  ❌ [Error] 유효한 숫자를 입력해주세요."
        fi
    done
else
    echo "❌ [Error] 프로젝트 폴더가 없습니다."
    exit 1
fi

# =====================================================================
# 2. 분석 카테고리 선택 (host, rRNA, viral 등)
# =====================================================================
STAR_DIR="$PROJECT_DIR/02-2.STAR_alignment"
if [ ! -d "$STAR_DIR" ]; then
    echo "❌ [Error] STAR alignment 디렉토리를 찾을 수 없습니다: $STAR_DIR"
    exit 1
fi

echo -e "\n[ 분석할 타겟 카테고리 선택 ]"
TARGET_DIRS=()
for d in "$STAR_DIR"/*/; do
    [ -d "$d" ] || continue
    TARGET_DIRS+=("$(basename "$d")")
done

if [ ${#TARGET_DIRS[@]} -gt 0 ]; then
    for i in "${!TARGET_DIRS[@]}"; do
        num=$((i+1))
        echo "  $num) ${TARGET_DIRS[$i]}"
    done
    echo ""
    while true; do
        read -p "▶ 카테고리 번호를 선택하세요: " cat_choice
        if [[ "$cat_choice" =~ ^[0-9]+$ ]] && [ "$cat_choice" -ge 1 ] && [ "$cat_choice" -le "${#TARGET_DIRS[@]}" ]; then
            SELECTED_CAT="${TARGET_DIRS[$((cat_choice-1))]}"
            TARGET_BAM_DIR="$STAR_DIR/$SELECTED_CAT"
            echo "  ▷ 선택된 카테고리: $SELECTED_CAT"
            break
        else
            echo "  ❌ [Error] 유효한 숫자를 입력해주세요."
        fi
    done
else
    echo "❌ [Error] $STAR_DIR 하위에 타겟 폴더가 없습니다."
    exit 1
fi

# 분석 결과 저장 디렉토리를 상위로 이동 배치 (권한 문제 방지)
ANALYSIS_DIR="$PROJECT_DIR/04_analysis/featureCounts/$SELECTED_CAT"
mkdir -p "$ANALYSIS_DIR"

# =====================================================================
# 3. BAM 파일 선택
# =====================================================================
echo -e "\n[ 분석할 BAM 파일 선택 ]"
BAM_FILES=()
for bam in "$TARGET_BAM_DIR"/*.bam; do
    [ -e "$bam" ] || continue
    BAM_FILES+=("$bam")
done

if [ ${#BAM_FILES[@]} -gt 0 ]; then
    for i in "${!BAM_FILES[@]}"; do
        num=$((i+1))
        echo "  $num) $(basename "${BAM_FILES[$i]}")"
    done
    echo ""
    while true; do
        read -p "▶ 분석할 BAM 파일 번호를 선택하세요: " bam_choice
        if [[ "$bam_choice" =~ ^[0-9]+$ ]] && [ "$bam_choice" -ge 1 ] && [ "$bam_choice" -le "${#BAM_FILES[@]}" ]; then
            TARGET_BAM="${BAM_FILES[$((bam_choice-1))]}"
            BAM_BASENAME="$(basename "$TARGET_BAM" .bam)"
            echo "  ▷ 선택된 BAM 파일: $BAM_BASENAME"
            break
        else
            echo "  ❌ [Error] 유효한 숫자를 입력해주세요."
        fi
    done
else
    echo "❌ [Error] $TARGET_BAM_DIR 에 BAM 파일이 존재하지 않습니다."
    exit 1
fi

# =====================================================================
# 4. GTF 파일 자동 스캔 및 선택 + 정제된 filtered.gtf 생성
# =====================================================================
echo -e "\n[ Reference GTF 파일 선택 ]"
GTF_FILES=($(find "$BASE_DIR/reference" -type f -name "*.gtf"))

if [ ${#GTF_FILES[@]} -gt 0 ]; then
    for i in "${!GTF_FILES[@]}"; do
        num=$((i+1))
        echo "  $num) ${GTF_FILES[$i]}"
    done
    echo ""
    while true; do
        read -p "▶ 사용할 GTF 파일 번호를 선택하세요: " gtf_choice
        if [[ "$gtf_choice" =~ ^[0-9]+$ ]] && [ "$gtf_choice" -ge 1 ] && [ "$gtf_choice" -le "${#GTF_FILES[@]}" ]; then
            TARGET_GTF="${GTF_FILES[$((gtf_choice-1))]}"
            echo "  ▷ 선택된 GTF 파일: $TARGET_GTF"
            break
        else
            echo "  ❌ [Error] 유효한 숫자를 입력해주세요."
        fi
    done
else
    echo "❌ [Error] GTF 파일이 존재하지 않습니다."
    exit 1
fi

# 정제된 GTF 파일 경로 설정 및 생성 (공백 gene_id 및 unknown_transcript 사전 차단)
FILTERED_GTF="$ANALYSIS_DIR/filtered_annotation.gtf"
echo "▶ GTF 파일 정제 중 (공백 gene_id 및 unknown_transcript 제거)..."
grep -v -E 'gene_id ""|unknown_transcript' "$TARGET_GTF" > "$FILTERED_GTF"
echo "✓ 정제된 GTF 생성 완료: $FILTERED_GTF"

# =====================================================================
# 5. 분석 기본 설정
# =====================================================================
echo -e "\n[ 분석 기본 설정 ]"
read -p "▶ Paired-end 데이터입니까? (y/n): " IS_PE
read -p "▶ 사용할 스레드 수를 입력하세요 (기본값: 8): " THREADS
THREADS=${THREADS:-8}
read -p "▶ 다중 매핑(Multimapping) 리드를 카운트에 포함하시겠습니까? (y/n): " COUNT_MULTIMAP

RAW_OUTPUT="$ANALYSIS_DIR/${BAM_BASENAME}_raw_counts.txt"
CLEAN_OUTPUT="$ANALYSIS_DIR/${BAM_BASENAME}_filtered_sorted_genes.txt"
DIST_OUTPUT="$ANALYSIS_DIR/${BAM_BASENAME}_category_distribution.txt"
FC_LOG="$ANALYSIS_DIR/${BAM_BASENAME}_featureCounts.log"

COMMON_FEATURES=("exon" "transcript" "gene" "직접 입력")
COMMON_ATTRIBUTES=("gene_id" "locus_tag" "직접 입력")
EXTRA_ATTRIBUTES=("gbkey" "gene_biotype" "사용 안 함 (추가 정보 없음)" "직접 입력")

# =====================================================================
# 6. 옵션 선택 및 featureCounts 실행
# =====================================================================
while true; do
    echo -e "\n[ featureCounts 세부 옵션 설정 ]"
    
    # Feature Type 설정
    echo "사용 가능한 Feature Type:"
    for i in "${!COMMON_FEATURES[@]}"; do echo "  $((i+1))) ${COMMON_FEATURES[$i]}"; done
    while true; do
        read -p "▶ Feature Type 번호를 선택하세요: " f_choice
        if [[ "$f_choice" =~ ^[0-9]+$ ]] && [ "$f_choice" -ge 1 ] && [ "$f_choice" -le "${#COMMON_FEATURES[@]}" ]; then
            if [ "$f_choice" -eq "${#COMMON_FEATURES[@]}" ]; then read -p "  ▷ 직접 입력: " FEATURE_TYPE
            else FEATURE_TYPE="${COMMON_FEATURES[$((f_choice-1))]}"; fi
            break
        else echo "  ❌ [Error] 유효한 숫자를 입력해주세요."; fi
    done

    # Primary Attribute Type 설정 (-g)
    echo -e "\n사용 가능한 기본 Attribute Type (그룹화 기준):"
    for i in "${!COMMON_ATTRIBUTES[@]}"; do echo "  $((i+1))) ${COMMON_ATTRIBUTES[$i]}"; done
    while true; do
        read -p "▶ Attribute Type 번호를 선택하세요: " a_choice
        if [[ "$a_choice" =~ ^[0-9]+$ ]] && [ "$a_choice" -ge 1 ] && [ "$a_choice" -le "${#COMMON_ATTRIBUTES[@]}" ]; then
            if [ "$a_choice" -eq "${#COMMON_ATTRIBUTES[@]}" ]; then read -p "  ▷ 직접 입력: " ATTRIBUTE_TYPE
            else ATTRIBUTE_TYPE="${COMMON_ATTRIBUTES[$((a_choice-1))]}"; fi
            break
        else echo "  ❌ [Error] 유효한 숫자를 입력해주세요."; fi
    done

    # Extra Attribute Type 설정 (--extraAttributes)
    echo -e "\n사용 가능한 추가 정보 Attribute Type (분포도 확인용):"
    for i in "${!EXTRA_ATTRIBUTES[@]}"; do echo "  $((i+1))) ${EXTRA_ATTRIBUTES[$i]}"; done
    while true; do
        read -p "▶ 추가 정보 번호를 선택하세요: " x_choice
        if [[ "$x_choice" =~ ^[0-9]+$ ]] && [ "$x_choice" -ge 1 ] && [ "$x_choice" -le "${#EXTRA_ATTRIBUTES[@]}" ]; then
            if [ "$x_choice" -eq $((${#EXTRA_ATTRIBUTES[@]} - 1)) ]; then EXTRA_ATTR=""
            elif [ "$x_choice" -eq "${#EXTRA_ATTRIBUTES[@]}" ]; then read -p "  ▷ 직접 입력: " EXTRA_ATTR
            else EXTRA_ATTR="${EXTRA_ATTRIBUTES[$((x_choice-1))]}"; fi
            break
        else echo "  ❌ [Error] 유효한 숫자를 입력해주세요."; fi
    done

    FC_OPTIONS="-T $THREADS -t $FEATURE_TYPE -g $ATTRIBUTE_TYPE"
    if [[ "$IS_PE" =~ ^[Yy]$ ]]; then FC_OPTIONS="$FC_OPTIONS -p"; fi
    if [[ "$COUNT_MULTIMAP" =~ ^[Yy]$ ]]; then FC_OPTIONS="$FC_OPTIONS -M"; fi
    if [ -n "$EXTRA_ATTR" ]; then FC_OPTIONS="$FC_OPTIONS --extraAttributes $EXTRA_ATTR"; fi

    echo -e "\n========================================================"
    echo "▶ featureCounts 실행 중..."
    echo "명령어: featureCounts $FC_OPTIONS -a $FILTERED_GTF -o $RAW_OUTPUT $TARGET_BAM"
    echo "========================================================"

    set +e
    featureCounts $FC_OPTIONS -a "$FILTERED_GTF" -o "$RAW_OUTPUT" "$TARGET_BAM" > "$FC_LOG" 2>&1
    FC_EXIT_CODE=$?
    set -e

    if [ $FC_EXIT_CODE -ne 0 ]; then
        echo -e "\n❌ [Error] featureCounts 실행 중 오류가 발생했습니다!"
        cat "$FC_LOG"
        echo "🔄 세부 옵션 설정 단계로 돌아갑니다."
        continue
    fi

    echo "✓ 원본 결과 저장 완료: $RAW_OUTPUT"
    break
done

# =====================================================================
# 7. 결과 정제 및 통계 요약 산출
# =====================================================================
echo -e "\n========================================================"
echo "▶ 결과 필터링 및 카테고리 통계 추출 중..."
echo "========================================================"

if [ -n "$EXTRA_ATTR" ]; then
    # 세미콜론으로 중복된 속성값(예: mRNA;mRNA;mRNA)을 고유값 하나로 정리하는 awk 처리
    tail -n +3 "$RAW_OUTPUT" | awk '$8 > 10 {
        n = split($7, arr, ";");
        delete seen;
        unique_attr = "";
        for(i=1; i<=n; i++) {
            if(!seen[arr[i]] && arr[i] != "") {
                seen[arr[i]] = 1;
                unique_attr = (unique_attr == "" ? arr[i] : unique_attr "," arr[i]);
            }
        }
        print $1 "\t" unique_attr "\t" $8
    }' | sort -k3,3nr > "$CLEAN_OUTPUT"

    # 카테고리 분포도 요약 (첫 번째 속성 기준 합산)
    tail -n +3 "$RAW_OUTPUT" | awk -F '\t' '{
        split($7, arr, ";");
        if(arr[1] != "") sum[arr[1]] += $8;
    } END {for (k in sum) print k "\t" sum[k]}' | sort -k2,2nr > "$DIST_OUTPUT"
    
    echo "✓ 개별 유전자 리스트 저장 완료: $(basename "$CLEAN_OUTPUT")"
    echo "✓ 카테고리별 분포 리스트 저장 완료: $(basename "$DIST_OUTPUT")"
else
    tail -n +3 "$RAW_OUTPUT" | awk '$7 > 10 {print $1 "\t" $7}' | sort -k2,2nr > "$CLEAN_OUTPUT"
    echo "✓ 개별 유전자 리스트 저장 완료: $(basename "$CLEAN_OUTPUT")"
    echo "⚠️ 추가 정보를 선택하지 않아 카테고리 분포도(Summary)는 생성되지 않았습니다."
fi

echo -e "\n🎉 모든 분석이 완료되었습니다."
