import os
import subprocess

# 构造参数列表
args = [
    "prepsubband_cu",  # 程序名 (相当于 argv[0])
    "-IOlog",          # 参数
    "-nobary",         # 参数
    "-numout",         # 参数
    "32768",           # 参数值
    "-nsub",           # 参数
    "3280",            # 参数值
    "-lodm",           # 参数
    "1000",            # 参数值
    "-dmstep",         # 参数
    "2.0",             # 参数值
    "-numdms",         # 参数
    "1000",            # 参数值
    "-downsamp",       # 参数
    "8",               # 参数值
    "-mask",           # 参数
    "/presto_data/test/FRB/rfi/test1_rfifind.mask",  # 参数值
    "-o",              # 参数
    "/presto_data/test/FRB/psb_tmp/FRB121102_tracking-M01_0706_ds1_0",  # 输出路径
    "/presto_data/test/FRB/FRB121102_tracking-M01_0706_ds1_0.fits"  # 输入文件
]

# 设置工作目录
work_dir = "/presto_data/test/FRB"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

# 循环执行 10 次
total_runs = 10
successful_runs = 0
failed_runs = 0

print(f"开始循环执行 {total_runs} 次 call_prepsubband_cu...")
print("=" * 60)

for run_num in range(1, total_runs + 1):
    print(f"\n--- 第 {run_num} 次执行 ---")
    # 删除结果目录（如果存在），确保干净环境
    result_dir = "/presto_data/test/FRB/psb_tmp"
    if os.path.exists(result_dir):
        try:
            import shutil
            shutil.rmtree(result_dir)
            print(f"已删除结果目录: {result_dir}")
        except Exception as e:
            print(f"删除结果目录 {result_dir} 时出错: {e}")
    else:
        print(f"结果目录 {result_dir} 不存在，无需删除")
    os.makedirs(result_dir, exist_ok=True)

    # 在子进程中调用动态库
    try:
        result = subprocess.run(
            ["python3", "-c", f"""
import ctypes
from ctypes import c_char_p, c_int

# 加载动态库
lib = ctypes.CDLL('/home/soft/presto/lib/libprepsubband_cu_lib.so')
lib.call_prepsubband_cu.argtypes = [c_int, ctypes.POINTER(c_char_p)]
lib.call_prepsubband_cu.restype = c_int

# 构造参数并调用
argc = {len(args)}
argv = (c_char_p * argc)(*map(lambda x: x.encode('utf-8'), {args}))
result = lib.call_prepsubband_cu(argc, argv)
if result != 0:
    raise RuntimeError(f'call_prepsubband_cu failed with return code {{result}}')
"""],
            cwd=work_dir,
            env=os.environ,
            capture_output=True,
            text=True,
        )
        print(result.stdout)
        print("✓ 执行成功!")
        successful_runs += 1
    except subprocess.CalledProcessError as e:
        print(f"✗ 执行失败，返回码: {e.returncode}")
        print(e.stderr)
        failed_runs += 1

# 输出总结
print("\n" + "=" * 60)
print("执行总结:")
print(f"总执行次数: {total_runs}")
print(f"成功次数: {successful_runs}")
print(f"失败次数: {failed_runs}")
print(f"成功率: {(successful_runs/total_runs)*100:.1f}%")