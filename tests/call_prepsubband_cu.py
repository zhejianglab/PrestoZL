import ctypes
from ctypes import CDLL, POINTER, c_char_p, c_int

# 加载共享库
lib = CDLL('/home/soft/presto/lib/libprepsubband_cu_lib.so')

# 定义 call_prepsubband_cu 的函数签名
lib.call_prepsubband_cu.argtypes = [c_int, POINTER(c_char_p)]
lib.call_prepsubband_cu.restype = c_int

# 构造参数列表，拆分为与命令行一致的形式
args = [
    "prepsubband_cu",  # 程序名 (相当于 argv[0])
    "-check",
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
    "/presto_data/test/FRB/test1_rfifind.mask",  # 参数值
    "-o",              # 参数
    "/presto_data/test/FRB/psb_tmp/FRB121102_tracking-M01_0706_ds1_0",  # 输出路径
    "/presto_data/test/FRB/FRB121102_tracking-M01_0706_ds1_0.fits"  # 输入文件
]

# 转换为 C 风格的参数
argc = len(args)
argv = (c_char_p * argc)(*map(lambda x: x.encode('utf-8'), args))

# 调用 call_prepsubband_cu
result = lib.call_prepsubband_cu(argc, argv)

# 输出结果
print(f"call_prepsubband_cu result: {result}")