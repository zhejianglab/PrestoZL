import ctypes
from ctypes import CDLL, POINTER, c_char_p, c_int

# 加载共享库
lib = CDLL('/home/soft/presto/lib/librealfft_lib.so')

# 定义 call_realfft 的函数签名
lib.call_realfft.argtypes = [c_int, POINTER(c_char_p)]
lib.call_realfft.restype = c_int

# 构造参数列表，拆分为与命令行一致的形式
args = [
    "realfft",  # 程序名 (相当于 argv[0])
    "-IOlog",     # 参数
    "-fwd",     # 参数
    "-outdir",  # 参数
    "/presto_data/test/FRB/psb_tmp",    # 参数值
    "/presto_data/test/FRB/psb/FRB121102_tracking-M01_0706_ds1_0_DM2998.00.dat"  # 输入文件
]

# 转换为 C 风格的参数
argc = len(args)
argv = (c_char_p * argc)(*map(lambda x: x.encode('utf-8'), args))

# 调用 call_realfft
result = lib.call_realfft(argc, argv)

# 输出结果
print(f"call_realfft result: {result}")