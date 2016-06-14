BEGIN{
	C["A"] = "T"; C["C"] = "G"; C["G"] = "C"; C["T"] = "A"
	C["a"] = "t"; C["c"] = "g"; C["g"] = "c"; C["t"] = "a"
}
(COMP > 0) {
	for (F = 1; F <= NF; F++) {
		SEQ=""
		for (i = length($F); i > 0; i--) {
			if (substr($F, i, 1) in C) {
				SEQ = SEQ C[substr($F, i, 1)]
			} else {
				SEQ = SEQ "N"
			}
		}
		$F = SEQ
	}
} (COMP == 0) {
	for (F = 1; F <= NF; F++) {
		SEQ=""
		for (i = length($F); i > 0; i--) {
			SEQ = SEQ substr($F, i, 1)
		}
		$F = SEQ
	}
}{
	print $0
}
