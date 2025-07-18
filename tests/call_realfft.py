import ctypes
import os
from ctypes import CDLL, POINTER, c_char_p, c_int

# 加载共享库
lib = CDLL('/home/soft/presto/lib/librealfft_lib.so')

# 定义 call_realfft 的函数签名
lib.call_realfft.argtypes = [c_int, POINTER(c_char_p)]
lib.call_realfft.restype = c_int

# 设置PRESTO环境变量
if "PRESTO" not in os.environ:
    os.environ["PRESTO"] = "/home/soft/presto"

# 切换到工作目录
working_dir = "/presto_data/test/FRB"
original_cwd = os.getcwd()
os.chdir(working_dir)

print(f"当前工作目录: {os.getcwd()}")

# 构造参数列表，与日志中失败的命令完全一致
args = [
    "realfft",  # 程序名 (相当于 argv[0])
    "-IOlog",     # 参数
    "-fwd",     # 参数
    "-outdir",  # 参数
    "/presto_data/test/FRB/psb_tmp",    # 参数值
    "/presto_data/test/FRB/psb/FRB121102_tracking-M01_0706_ds1_0_DM2998.00.dat"  # 输入文件
]

print(f"执行命令: {' '.join(args)}")

# 检查输入文件是否存在
input_file = args[-1]
if not os.path.exists(input_file):
    print(f"错误: 输入文件不存在: {input_file}")
    exit(1)

# 检查输出目录是否存在
output_dir = args[-2]
if not os.path.exists(output_dir):
    print(f"错误: 输出目录不存在: {output_dir}")
    exit(1)

print(f"输入文件存在: {input_file}")
print(f"输出目录存在: {output_dir}")

# 转换为 C 风格的参数
argc = len(args)
argv = (c_char_p * argc)(*map(lambda x: x.encode('utf-8'), args))

try:
    # 调用 call_realfft
    print("开始执行 realfft...")
    result = lib.call_realfft(argc, argv)
    
    # 输出结果
    print(f"call_realfft result: {result}")
    
    if result == 0:
        print("realfft 执行成功!")
    else:
        print(f"realfft 执行失败，返回码: {result}")
        
except Exception as e:
    print(f"执行过程中发生异常: {e}")
finally:
    # 恢复原始工作目录
    os.chdir(original_cwd)
    print(f"恢复工作目录: {os.getcwd()}") 