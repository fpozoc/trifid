.PHONY: build check clean-pyc clean-test clean docker-bash docker-build docker-info docker-stats env help lint site sort test type
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate
PACKAGE_NAME = trifid

build:
	@echo Testing and then building
	python setup.py sdist bdist_wheel

check: format lint sort test type
clean-pyc:
	@echo Cleaning up python cache files
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	@echo Cleaning test files
	rm -f .coverage
	rm -f .coverage.*
	find . -name '.pytest_cache' -exec rm -fr {} +

clean: clean-pyc clean-test
	@echo Cleaning build artifacts
	find . -name '.my_cache' -exec rm -fr {} +
	rm -rf logs/

docker-bash:
	@echo Running bash in the docker container
	docker exec -it trifid bash

docker-build:
	@echo Building the docker container
	docker build -t trifid .

docker-info:
	@echo Displaying docker info
	docker info

docker-stats:
	@echo Displaying docker stats, images, containers and volumes
	docker images -a
	@echo " "
	docker ps -a $(DOCKER_PS_ARGS)
	@echo " "
	docker volume ls

env:
	@echo Creating the conda environment
	git init
	mamba env create -f environment.yml
	($(CONDA_ACTIVATE) $(PACKAGE_NAME))
	mamba run -n $(PACKAGE_NAME) pre-commit install

format:
	@echo Formatting code
	black src/$(PACKAGE_NAME)

help:
	@echo Show help on available commands
	grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

lint:
	@echo Linting with pylint
	mamba run -n $(PACKAGE_NAME) pylint src/$(PACKAGE_NAME) --reports=y

site:
	@echo Building the documentation
	mamba run -n $(PACKAGE_NAME) mike deploy --update-aliases 0.0 latest
	mamba run -n $(PACKAGE_NAME) mike set-default latest

sort:
	@echo Sorting imports with isort
	mamba run -n $(PACKAGE_NAME) isort src/$(PACKAGE_NAME)

test:
	@echo Running tests
	@echo Testing
	pytest --cov=src/$(PACKAGE_NAME) tests/

type:
	@echo Running mypy
	mamba run -n $(PACKAGE_NAME) mypy src/$(PACKAGE_NAME)
