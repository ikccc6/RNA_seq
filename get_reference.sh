#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <Accession_ID>"
    exit 1
fi

ACCESSION=$1
BASE_DIR="$HOME/reference/genome"
OUT_DIR="$BASE_DIR/$ACCESSION"

mkdir -p "$OUT_DIR"

if [[ "$ACCESSION" == GCF_* ]] || [[ "$ACCESSION" == GCA_* ]]; then
    # 1. 기존 파일 존재 여부 확인
    NEED_FASTA=true
    NEED_GTF=true
    
    if [ -f "$OUT_DIR/${ACCESSION}.fna" ]; then
        echo "[Info] FASTA file already exists."
        NEED_FASTA=false
    fi
    
    if [ -f "$OUT_DIR/${ACCESSION}.gtf" ]; then
        echo "[Info] GTF file already exists."
        NEED_GTF=false
    fi

    # 2. 둘 다 존재하면 스크립트 종료
    if [ "$NEED_FASTA" = false ] && [ "$NEED_GTF" = false ]; then
        echo "[Info] All files exist. Skipping download."
        exit 0
    fi

    # 3. 필요한 파일에 맞게 API URL 구성
    API_URL="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${ACCESSION}/download?"
    
    if [ "$NEED_FASTA" = true ]; then
        API_URL="${API_URL}include_annotation_type=GENOME_FASTA&"
    fi
    if [ "$NEED_GTF" = true ]; then
        API_URL="${API_URL}include_annotation_type=GENOME_GTF"
    fi
    
    # URL 끝에 '&'가 남아있다면 제거
    API_URL=${API_URL%&}

    echo "[Info] Downloading missing files..."
    ZIP_FILE="$OUT_DIR/${ACCESSION}.zip"
    
    wget -q --show-progress -O "$ZIP_FILE" "$API_URL"
    
    if [ $? -eq 0 ]; then
        echo "[Info] Extracting files..."
        unzip -q "$ZIP_FILE" -d "$OUT_DIR"
        
        # .fna 이동 (필요했던 경우에만)
        if [ "$NEED_FASTA" = true ] && ls "$OUT_DIR/ncbi_dataset/data/$ACCESSION/"*.fna 1> /dev/null 2>&1; then
            mv "$OUT_DIR/ncbi_dataset/data/$ACCESSION/"*.fna "$OUT_DIR/${ACCESSION}.fna"
        fi
        
        # .gtf 이동 (필요했던 경우에만)
        if [ "$NEED_GTF" = true ]; then
            if ls "$OUT_DIR/ncbi_dataset/data/$ACCESSION/"*.gtf 1> /dev/null 2>&1; then
                mv "$OUT_DIR/ncbi_dataset/data/$ACCESSION/"*.gtf "$OUT_DIR/${ACCESSION}.gtf"
            else
                echo "[Warning] No GTF file found for this assembly."
            fi
        fi
        
        # 정리
        rm -rf "$OUT_DIR/ncbi_dataset" "$OUT_DIR/README.md" "$ZIP_FILE"
        echo "[Success] Done: $OUT_DIR"
    else
        echo "[Error] Download failed."
        rm -rf "$ZIP_FILE"
        exit 1
    fi

else
    # 단일 시퀀스 로직
    if [ -f "$OUT_DIR/${ACCESSION}.fna" ]; then
        echo "[Info] FASTA file already exists. Skipping download."
        exit 0
    fi

    echo "[Info] Single Sequence Accession detected. Downloading FASTA..."
    wget -q --show-progress -O "$OUT_DIR/${ACCESSION}.fna" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACCESSION}&rettype=fasta&retmode=text"
    
    if [ -s "$OUT_DIR/${ACCESSION}.fna" ]; then
        echo "[Success] Done: $OUT_DIR/${ACCESSION}.fna"
    else
        echo "[Error] Download failed or invalid Accession ID."
        rm -rf "$OUT_DIR"
        exit 1
    fi
fi
