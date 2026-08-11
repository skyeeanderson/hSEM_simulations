#!/bin/bash
SCENARIOS=(
  abundant_50   abundant_500   abundant_1000
  rare_50       rare_500       rare_1000
  abundant_50_meso   abundant_500_meso   abundant_1000_meso
  rare_50_meso       rare_500_meso       rare_1000_meso
)
LANDSCAPES=(TRUE FALSE)
mkdir -p /scratch/user/uqsand24/simulations/slurm_v4rn

get_time() {
  case $1 in
    *_50|*_50_meso)     echo "6:00:00"  ;;
    *_500|*_500_meso)   echo "16:00:00" ;;
    *_1000|*_1000_meso) echo "36:00:00" ;;
  esac
}
get_mem() {
  case $1 in
    *_50|*_50_meso)     echo "8GB"  ;;
    *_500|*_500_meso)   echo "16GB" ;;
    *_1000|*_1000_meso) echo "32GB" ;;
  esac
}

for SC in "${SCENARIOS[@]}"; do
  TIME=$(get_time $SC)
  MEM=$(get_mem $SC)
  for LAND in "${LANDSCAPES[@]}"; do
    TAG=$([ "$LAND" == "TRUE" ] && echo "LT" || echo "LF")
    OUT="/scratch/user/uqsand24/simulations/slurm_v4rn/run_${SC}_${TAG}.sh"
    cat > "$OUT" << SLURM
#!/bin/bash
#SBATCH --job-name=v4rn_${SC}_${TAG}
#SBATCH --output=/scratch/user/uqsand24/simulations/logs/v4rn_${SC}_${TAG}_%a.out
#SBATCH --error=/scratch/user/uqsand24/simulations/logs/v4rn_${SC}_${TAG}_%a.err
#SBATCH --array=1-100%10
#SBATCH --time=${TIME}
#SBATCH --mem=${MEM}
#SBATCH --cpus-per-task=1
#SBATCH --account=a_senv_bcs

module load rjags/4-17-foss-2024a-r-4.4.2

export SETTING=HPC
export SCENARIO=${SC}
export LANDSCAPE=${LAND}

mkdir -p /scratch/user/uqsand24/simulations/results/v4rn_${SC}_${TAG}
mkdir -p /scratch/user/uqsand24/simulations/logs

Rscript /scratch/user/uqsand24/simulations/scripts/01_simulate_and_fit.R
SLURM
    chmod +x "$OUT"
  done
done
echo "Generated $(ls /scratch/user/uqsand24/simulations/slurm_v4rn/*.sh | wc -l) scripts"
