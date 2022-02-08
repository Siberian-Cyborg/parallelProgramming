def read_all_text(p):
    with open(p, "r") as f:
        return f.read()

def process_submission(src, dest, headers, replace_search=None, replace=None):
    src_txt = read_all_text(src)
    for h in headers:
        header_txt = f"// ---------------- BEGIN: {h} ----------------\n{read_all_text(h)}\n// ----------------  END:  {h} ----------------\n"
        guard = f'#include "{h}"'
        #src_txt = "#define SUBMISSION\n\n" + src_txt
        src_txt = src_txt.replace(guard, header_txt)
        if replace_search is not None and replace is not None:
            src_txt = src_txt.replace(replace_search, replace)
    with open(dest, 'w') as f:
        f.write(src_txt)

if __name__ == "__main__":
    process_submission("hybridbh.cpp", "submission/hybridbh_submission.txt", ["hybridbh.h"])
    process_submission("bonusge.cpp", "submission/bonusge_submission.txt", ["bonusge.h"])
    #for t in range(1, 33):
    #    for s in ['static', 'guided']:
    #        process_submission("ompbh.cpp", f"submission/omp_submission_{t}_{s[0]}.cpp", ["ompbh.h"], "schedule(guided) num_threads(24)", f"schedule({s}) num_threads({t})")
