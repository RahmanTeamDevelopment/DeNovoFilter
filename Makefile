unittest:
	echo "Unit tests:"
	source env/bin/activate; \
	pytest test/unit -v
	
smoketest:
	test/smoke/test_installation.sh
	
test: smoketest unittest

