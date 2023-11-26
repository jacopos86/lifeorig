install :
	python setup.py install
.PHONY :
	clean test
clean :
	rm -rf ./lifeorig/*~ ./lifeorig/__pycache__ ./build/lib/lifeorig/* ./__pycache__
