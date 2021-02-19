for f in *.xls; do
	mv -- "$f" "${f%.xls}C.csv"
done
