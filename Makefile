# Makefile for qutip-bo

redb = "\033[1;31m"$1"\033[0m"
greenb = "\033[1;32m"$1"\033[0m"
green = "\033[0;32m"$1"\033[0m"
yellow = "\033[0;33m"$1"\033[0m"
whiteb = "\033[1;37m"$1"\033[0m"
DOCKERHUB_ORGANISATION=julianiacoponi
VERSION?=$(shell git branch | grep "*" | sed -e "s/* //")

help:
	@echo $(VERSION)
	@echo $(call greenb,"Please use make <target>' where <target> is one of:")
	@echo $(call redb," NOTHING IMPLEMENTED PLEASE DON'T USE!")" (please...)"

pull-%:
	docker pull $(DOCKERHUB_ORGANISATION)/$*
export_version:
	@echo "__version__ = '$(VERSION)'" > qutip/__init__.py
	@echo $(VERSION) > .VERSION

qutip: export_version
	docker build --tag=${DOCKERHUB_ORGANISATION}/$@:$(VERSION) -f dockerfiles/Dockerfile.qutip-bo .
	docker tag ${DOCKERHUB_ORGANISATION}/$@:${VERSION} $@
	docker tag ${DOCKERHUB_ORGANISATION}/$@:${VERSION} ${DOCKERHUB_ORGANISATION}/$@

push-%:
	docker tag $(DOCKERHUB_ORGANISATION)/$*:$(VERSION) $(DOCKERHUB_ORGANISATION)/$*:latest
	docker push $(DOCKERHUB_ORGANISATION)/$*:$(VERSION)
	docker push $(DOCKERHUB_ORGANISATION)/$*:latest

qutip-dev: export_version
	@echo  $(call green,"Launch qutip-dev image")
	docker build --tag=qutip-dev -f Dockerfile.qutip-dev .
	@# touch the file so that docker doesn't make a directory for it
	@touch $(shell pwd)/.dev-image-history
	docker run -it -v $(shell pwd):/qcbo/ -v $(shell pwd)/.dev-image-history:/root/.bash_history --entrypoint=/bin/bash qutip-dev

gpyopt-dev: export_version
	@echo  $(call green,"Launch gpyopt-dev image")
	docker build --tag=gpyopt-dev -f Dockerfile.gpyopt-dev .
	@# touch the file so that docker doesn't make a directory for it
	@touch $(shell pwd)/.dev-image-history
	docker run -it -v $(shell pwd):/qcbo/ -v $(shell pwd)/.dev-image-history:/root/.bash_history --entrypoint=/bin/bash gpyopt-dev

qcbo-dev: export_version
	@echo  $(call green,"Launch qcbo-dev image")
	docker build --tag=qcbo-dev -f Dockerfile.qcbo-dev .
	@# touch the file so that docker doesn't make a directory for it
	@touch $(shell pwd)/.dev-image-history
	docker run -it -v $(shell pwd):/qcbo/ -v $(shell pwd)/.dev-image-history:/root/.bash_history --entrypoint=/bin/bash qcbo-dev
