import os
import subprocess

# 构造参数列表
args = [
    "prepfold",  # 程序名 (相当于 argv[0])
    "-IOlog",          # 参数
    "-cache",         # 参数
    "-json",         # 参数
    "-fine",           # 参数值
    "-topo",           # 参数
    "-npart",            # 参数值
    "256",           # 参数
    "-nsub",            # 参数值
    "164",         # 参数
    "-noxwin",             # 参数值
    "-accelcand",         # 参数
    "7",            # 参数值
    "-accelfile",       # 参数
    "./psb/FRB121102_tracking-M01_0706_ds1_0_DM187.00_ACCEL_50_JERK_60.cand",               # 参数值
    "-mask",           # 参数
    "./rfi/test1_rfifind.mask",  # 参数值
    "-p",              # 参数
    "0.000999876",  # 参数
    "-dm",  # 参数
    "187.00",  # 参数
    "-o",  # 参数
    "./pic_zj_cache/FRB121102_tracking-M01_0706_ds1_0_11",  # 参数
    "./FRB121102_tracking-M01_0706_ds1_0.fits"  # 输入文件
]

# 设置工作目录
work_dir = "/presto_data/test/FRB"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

# 循环执行 10 次
total_runs = 10
successful_runs = 0
failed_runs = 0

print(f"开始循环执行 {total_runs} 次 call_prepfold...")
print("=" * 60)

for run_num in range(1, total_runs + 1):
    print(f"\n--- 第 {run_num} 次执行 ---")
    # 删除结果目录（如果存在），确保干净环境
    result_dir = "/presto_data/test/FRB/pic_zj_cache"
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
lib = ctypes.CDLL('/home/soft/presto/lib/libprepfold_lib.so')
lib.call_prepfold.argtypes = [c_int, ctypes.POINTER(c_char_p)]
lib.call_prepfold.restype = c_int

# 构造参数并调用
argc = {len(args)}
argv = (c_char_p * argc)(*map(lambda x: x.encode('utf-8'), {args}))
result = lib.call_prepfold(argc, argv)
if result != 0:
    raise RuntimeError(f'call_prepfold failed with return code {{result}}')
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