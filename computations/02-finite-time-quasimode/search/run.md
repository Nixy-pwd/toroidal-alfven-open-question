# 复算候选搜索

Python 3.12、NumPy、SciPy。以下命令均在此目录中运行。高阶搜索使用较多内存。独立的有理数核验在上一级目录运行，只需标准库。

    export OPENBLAS_NUM_THREADS=1
    ALF_DROP=0 python continue_modes.py 24 40 56 --out run
    ALF_DROP=0 python continue_modes.py 72 88 --tau 1e-12 1e-13 2e-14 --out run
    python continue_modes.py 112 144 --tau 1e-12 1e-13 2e-14 --out run
    OPENBLAS_NUM_THREADS=2 python continue_modes.py 176 --tau 4e-14 --out run
    python enrich.py run 144 176
    python budget.py run 24 40 56 72 88 112 144 176
    python export.py 176 2490 run
    python ../consequences.py run/trial_176_2490.json

continue_modes.py 用低阶候选作为固定初始方向，解三次带梯度惩罚的反迭代以生成试探向量。预算优化允许在这些向量张成的空间中调整组合，只报告可达到的残差。p=176 空间含一个独立生成的高阶试探和三个嵌入的 p=144 试探。search/raw-search.csv 各列依次为阶数、梯度惩罚系数、物理空间的残差上界估计、实测三维梯度能量、与初始方向的近似质量内积绝对值。表中的积分是浮点求积。独立的 trial.json 把自由系数四舍五入到分母为 2 的 72 次方的有理数，再精确组装边界相切的通量。它的严格上界见上一级目录的核验输出。
