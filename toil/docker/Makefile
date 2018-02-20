# Definitions
build_output = runtime/margin_phase
runtime_fullpath = $(realpath runtime)
build_tool = runtime-container.DONE
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
name = quay.io/ucsc_cgl/margin_phase
tag = ${git_commit}

# Steps
build:
	# so that our docker tag and build match the current git hash and repo
	-rm -rf ../../../tempMarginPhase
	-rm -rf tempMarginPhase
	cp -r ../../../marginPhase ../../../tempMarginPhase
	-rm -rf ../../../tempMarginPhase/toil
	-rm -rf ../../../tempMarginPhase/cmake-build-debug
	-rm -rf ../../../tempMarginPhase/CMakeCache.txt
	-rm -rf ../../../tempMarginPhase/CMakeFiles
	-rm -rf ../../../tempMarginPhase/*.cmake
	-cd ../../../tempMarginPhase && make clean
	mv ../../../tempMarginPhase tempMarginPhase
	# do the docker build
	docker build -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest
	rm -rf tempMarginPhase
	touch ${build_tool}

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}
	docker push ${name}:latest

test: build
	python test.py

clean:
	-rm ${build_tool}
	-rm ${build_output}
