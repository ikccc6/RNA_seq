import os
import io
import subprocess
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def select_menu(items, prompt_text):
    for i, item in enumerate(items, 1):
        name = item.name if isinstance(item, Path) else item
        print(f"[{i}] {name}")

    while True:
        try:
            idx = int(input(prompt_text)) - 1
            if 0 <= idx < len(items):
                return items[idx]
            else:
                print(f"[Error] : 1에서 {len(items)} 사이의 번호를 입력해주세요.\n")
        except ValueError:
            print("[Error] : 문자가 아닌 숫자를 입력해주세요.\n")

# 1. Base 경로 설정
base_dir = Path(os.path.expanduser("~/RNA_seq"))
ref_base_dir = base_dir / "reference" / "genome"

# 2. Plot 모드 및 샘플 개수 입력 받기
print("Select Plot Mode:")
print("[1] Overlay")
print("[2] Separate rows")
mode_choice = input("Enter mode number: ").strip()
is_overlay = (mode_choice == '1')

num_bams = int(input("\nHow many samples (BAM files) to plot? "))
bams_info = []

for i in range(num_bams):
    print(f"\n=============================================")
    print(f"       Select Data for Sample {i+1}       ")
    print(f"=============================================")
    projects = sorted([p for p in base_dir.iterdir() if p.is_dir() and p.name[:6].isdigit()])
    sel_project = select_menu(projects, f"Select Project number for Sample {i+1}: ")
    
    align_dir = sel_project / "02-2.STAR_alignment"
    targets = sorted([d for d in align_dir.iterdir() if d.is_dir()])
    
    # [Step 1] 전체 Read 수 계산을 위한 BAM 선택
    print("\n[Step 1] Total Reads(분모) 계산용 BAM 파일 선택")
    sel_tot_target = select_menu(targets, "Select Target directory number: ")
    tot_bams = sorted(list(sel_tot_target.glob("*.bam")))
    sel_tot_bam = select_menu(tot_bams, "Select BAM file number: ")
    
    print("Calculating total input reads...")
    total_reads = int(subprocess.check_output(f"samtools view -c -F 256 {sel_tot_bam}", shell=True))
    print(f"-> Total reads calculated: {total_reads:,}")
    
    # [Step 2] 실제 그래프를 그릴 BAM 선택
    print("\n[Step 2] Coverage Depth(분자) 그래프용 BAM 파일 선택")
    sel_cov_target = select_menu(targets, "Select Target directory number: ")
    cov_bams = sorted(list(sel_cov_target.glob("*.bam")))
    sel_cov_bam = select_menu(cov_bams, "Select BAM file number: ")
    
    # Overlay 모드일 경우 범례에 표시할 이름 입력
    label_name = ""
    if is_overlay:
        label_name = input(f"\nEnter Legend Label for Sample {i+1}: ")
    
    bams_info.append({
        'project': sel_project,
        'bam': sel_cov_bam,
        'total_reads': total_reads,
        'label': label_name
    })

# 3. Reference Directory 선택 및 FASTA 파일 스캔
print("\n--- Select Reference ---")
ref_dirs = sorted([d for d in ref_base_dir.iterdir() if d.is_dir()])
sel_ref_dir = select_menu(ref_dirs, "Select Reference directory number: ")
fastas = sorted(list(sel_ref_dir.glob("*.fasta")) + list(sel_ref_dir.glob("*.fna")))

# 4. 파일 내부의 모든 Contig ID 추출 및 표시할 이름 지정
contig_ids = []
for fasta in fastas:
    with open(fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                contig_ids.append(line.strip().lstrip('>').split()[0])

ref_names = []
print("\n--- Name Reference Segments ---")
for cid in contig_ids:
    ref_name = input(f"Enter plot name for contig '{cid}' (e.g., PB2): ")
    ref_names.append(ref_name)

# 5. 데이터 추출 및 연산 (Contig ID 단위 순회)
data_matrix = []
for info in bams_info:
    row_data = []
    for cid in contig_ids:
        cmd = f"samtools depth -a -r {cid} {info['bam']}"
        depth_out = subprocess.check_output(cmd, shell=True).decode()
        
        if not depth_out.strip():
            row_data.append(pd.DataFrame(columns=['pos', 'norm_depth']))
            continue
            
        df = pd.read_csv(io.StringIO(depth_out), sep='\t', names=['contig', 'pos', 'depth'])
        df['norm_depth'] = df['depth'] / info['total_reads']
        row_data.append(df)
    data_matrix.append(row_data)

# 6. 다중 Subplots 생성
# Overlay 모드면 행(Row) 개수는 1개로 고정, 분리 모드면 샘플 개수만큼 할당
num_rows = 1 if is_overlay else num_bams
num_cols = len(contig_ids)

base_fontsize = 10 + (num_cols - 1) * 0.4
tick_fontsize = base_fontsize - 2

# ref_name(contig) 개수에 따른 유동적인 너비 계산
# 8개 기준 1.3 인치가 되도록 10.4를 개수로 나눔
# 개수가 적을 때 지나치게 넓어지는 것을 방지하기 위해 최대 2.5 인치로 제한
plot_width = min(2.5, 10.4 / num_cols)
plot_height = 2.5

fig, axes = plt.subplots(
    num_rows, 
    num_cols, 
    figsize=(plot_width * num_cols, plot_height * num_rows),
    constrained_layout=True)

# axes 배열 차원 일관성 유지 (1차원 혹은 단일 객체일 경우 2차원으로 변환)
if getattr(axes, 'ndim', 0) == 0:
    axes = np.array([[axes]])
elif axes.ndim == 1:
    axes = axes[np.newaxis, :] if num_rows == 1 else axes[:, np.newaxis]

# Y축, X축 최대 수치 연산
all_max_y = []
for b_idx in range(num_bams):
    for c in range(num_cols):
        if not data_matrix[b_idx][c].empty:
            all_max_y.append(data_matrix[b_idx][c]['norm_depth'].max())
global_max_y = max(all_max_y) if all_max_y else 0.01

def get_nice_max(val):
    if val == 0: return 0.025
    power = math.floor(math.log10(val))
    factor = 10 ** power
    ceil_val = math.ceil(val / factor) * factor
    return ceil_val if ceil_val > val else ceil_val + factor / 2

nice_max_y = get_nice_max(global_max_y)
y_ticks = [nice_max_y * i / 5 for i in range(6)]

col_max_x = []
for c in range(num_cols):
    max_x = max([data_matrix[b_idx][c]['pos'].max() for b_idx in range(num_bams) if not data_matrix[b_idx][c].empty], default=1000)
    col_max_x.append(max_x)

# 기본 색상 팔레트 설정 (Overlay 시 다채로운 색상 적용)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# 7. 그래프 렌더링
for r in range(num_rows):
    for c in range(num_cols):
        ax = axes[r, c]
        
        if is_overlay:
            # 여러 BAM 파일을 하나의 subplot(ax)에 겹쳐서 그림
            for b_idx in range(num_bams):
                df = data_matrix[b_idx][c]
                if not df.empty:
                    ax.plot(df['pos'], df['norm_depth'], color=colors[b_idx % len(colors)], linewidth=1.2, label=bams_info[b_idx]['label'])
        else:
            # 분리 모드일 경우 각 행(r)에 매칭되는 데이터를 단일 색(검정)으로 그림
            df = data_matrix[r][c]
            if not df.empty:
                ax.plot(df['pos'], df['norm_depth'], color='black', linewidth=1.2)
        
        # 첫 번째 행에만 Contig 이름 표시
        if r == 0:
            ax.set_title(ref_names[c], pad=15, fontsize=base_fontsize)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        ax.set_xlim(0, col_max_x[c])
        ax.margins(x=0, y=0)
        ax.set_ylim(0, nice_max_y)
        
        # 하단 그래프(또는 Overlay 단일 그래프)에만 X축 수치 표기
        if r == num_rows - 1:
            ax.set_xticks([0, col_max_x[c]])
        else:
            ax.set_xticks([])
            
        # 첫 번째 열에만 Y축 수치 표기
        if c == 0:
            ax.set_yticks(y_ticks)
        else:
            ax.set_yticks(y_ticks)
            ax.tick_params(labelleft=False)

        ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)
            
        # Overlay 모드이고 마지막 열(가장 우측)일 때만 외곽에 범례(Legend) 표시
        if is_overlay and c == num_cols - 1:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, prop={'size':6, 'family': 'sans-serif'}, handlelength=1.0)

# 투명 보조 축을 이용해 다중 행/단일 행 관계없이 Y축 길이를 기준으로 정중앙 라벨 정렬
gs = axes[0, 0].get_gridspec()
fig.supxlabel("\nNucleotide position (5' -> 3')", fontsize=base_fontsize)
fig.supylabel("Coverage Depth\n", fontsize=base_fontsize)

# 9. 최초 입력받은 Project의 디렉토리에 저장
out_dir = bams_info[0]['project'] / "05_visualization"
out_dir.mkdir(parents=True, exist_ok=True)

plot_type = "overlay" if is_overlay else "separate"
output_filepath = out_dir / f"Combined_{num_bams}_samples_{sel_ref_dir.name}_{plot_type}_coverage.png"
plt.savefig(output_filepath, dpi=300, bbox_inches='tight') # bbox_inches='tight'가 외부 범례 잘림 방지
plt.close()
print(f"\nSaved combined multi-plot: {output_filepath}")
