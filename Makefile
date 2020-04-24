init:
	pip install -r requirements.txt

test:
	py.base tests

.PHONY: init test
