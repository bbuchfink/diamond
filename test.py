# import sys
# sys.path.append("build")
# https://zhuanlan.zhihu.com/p/164598060?utm_id=0#:~:text=%E6%9E%84%E5%BB%BA%20Python%20C%20%E6%89%A9%E5%B1%95%E6%A8%A1%E5%9D%97%201%20%E5%9C%A8%20Python%20%E5%86%85%E9%83%A8%E6%89%A7%E8%A1%8C,%E5%9C%A8%20C%20%E4%B8%AD%E5%AE%9A%E4%B9%89%E5%85%A8%E5%B1%80%E5%B8%B8%E9%87%8F%EF%BC%8C%E5%B9%B6%E5%9C%A8%20Python%20%E4%B8%AD%E8%AE%BF%E9%97%AE%E5%AE%83%E4%BB%AC%205%20%E6%89%93%E5%8C%85%E5%92%8C%E5%8F%91%E5%B8%83Python%20C%E6%89%A9%E5%B1%95%E6%A8%A1%E5%9D%97
# https://messense.me/build-python-c-extensions-with-cmake#:~:text=%E4%BD%BF%E7%94%A8%20CMake%20%E6%9E%84%E5%BB%BA%20Python%20C%2FC%2B%2B%20%E6%89%A9%E5%B1%95%201%20%E8%B0%83%E7%A0%94%E5%8F%AF%E9%80%89%E6%96%B9%E6%A1%88,Cython%20%E4%B8%AD%E4%BD%BF%E7%94%A8%20abseil-cpp%20%E7%9A%84%20containers%20%E7%9A%84%E6%96%87%E7%AB%A0%EF%BC%8Cstay%20tuned.%20
from diamondpy.libdiamond import run, version
print(version())
try:
    run(
        "/home/hongning/.annopro/data/cafa4.dmnd", 
        "/home/hongning/AnnoPRO/test_proteins.fasta", 
        "blastp_output", 
        "blastp")
except RuntimeError as e:
    print(e)
print("done")