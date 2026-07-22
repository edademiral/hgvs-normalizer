from parse_coord import parse_variant

SV_MIN_BP = 50
CNV_MIN_BP = 1000

def variant_size(v):
    if v["size"]:
        num = v["size"].lower().replace("kb", "")
        try:
            return float(num) * 1000
        except ValueError:
            pass
    if v["start"] and v["end"]:
        return abs(int(v["end"]) - int(v["start"]))
    if v["ref"] and v["alt"]:
        return abs(len(v["ref"]) - len(v["alt"]))
    return None

def classify(v):
    if v["status"] == "needs_review":
        return "review"
    size = variant_size(v)
    if v["vtype"] == "SNV":
        return "SNV"
    if size is None:
        return "review"
    if size >= CNV_MIN_BP:
        return "CNV"
    elif size >= SV_MIN_BP:
        return "SV"
    else:
        return "SV"

if __name__ == "__main__":
    groups = {"SNV": [], "SV": [], "CNV": [], "review": []}
    with open("variants.tsv", encoding="utf-8") as f:
        for line in f:
            v = parse_variant(line)
            if v is None:
                continue
            grp = classify(v)
            groups[grp].append(v)
    for grp_name in ["SNV", "SV", "CNV", "review"]:
        items = groups[grp_name]
        print("=" * 50)
        print("  " + grp_name + " grubu  (" + str(len(items)) + " adet)")
        print("=" * 50)
        for v in items:
            size = variant_size(v)
            size_str = str(int(size)) + " bp" if size else "?"
            konum = str(v["chrom"]) + ":" + str(v["start"])
            if v["end"]:
                konum += "-" + str(v["end"])
            print("  " + konum + "  [" + str(v["vtype"]) + ", " + size_str + "]")
            print("     <- " + v["raw"].strip())
        print()
