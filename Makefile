fulltest:
	pytest -v -s --profile-svg --profile
	for i in prof/*.prof; do gprof2dot -f pstats $i | dot -Tsvg -o $i.svg; done

litetest:
	pytest  -m "not full" -v -s --profile-svg --profile
	for i in prof/*.prof; do gprof2dot -f pstats $i | dot -Tsvg -o $i.svg; done

singletest:
	cd ~/git_repos/_packages_/python/soxspipe
	time pytest -k "test_horne_extraction_function" -v -s --profile-svg --profile
	for i in prof/*.prof; do gprof2dot -f pstats $i | dot -Tsvg -o $i.svg; done

