;; Set up build environment using GNU Guix packages
;;
;; CC0 license, Pjotr Prins (c) 2022-2023
;;
;; To use this file to build HEAD:
;;
;;   guix build -f guix.scm
;;
;; To get a development container (emacs shell will work)
;;
;;   guix shell -C -D -f guix.scm
;;
;; For the tests you need /bin/bash. In a container create it with
;;
;;   mkdir -p /bin ; ln -s $GUIX_ENVIRONMENT/bin/bash /bin/bash
;;
;; To find tools
;;
;;   cd build
;;   cmake .. -DCMAKE_MAKE_PROGRAM=make -DCMAKE_C_COMPILER=gcc
;;   cmake --build . --verbose
;;
;; and run the tests with
;;
;;   env CC=gcc make
;;   ./tests/wfa.utest.sh

(use-modules
 ((guix licenses) #:prefix license:)
  (guix gexp)
  (guix packages)
  (guix git-download)
  (guix build-system cmake)
  (gnu packages algebra)
  (gnu packages autotools)
  (gnu packages base)
  (gnu packages bash)
  (gnu packages compression)
  (gnu packages build-tools)
  (gnu packages check)
  (gnu packages curl)
  (gnu packages gcc)
  (gnu packages gdb)
  (gnu packages llvm)
  (gnu packages parallel)
  (gnu packages pkg-config)
  (srfi srfi-1)
  (ice-9 popen)
  (ice-9 rdelim))

(define %source-dir (dirname (current-filename)))

(define %git-commit
    (read-string (open-pipe "git show HEAD | head -1 | cut -d ' ' -f 2" OPEN_READ)))

(define-public wfa2-lib-git
  (package
    (name "wfa2-lib-git")
    (version (git-version "1.3" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #f))
    (build-system cmake-build-system)
    (inputs
     `(("bash" ,bash)
       ("gdb" ,gdb)))
    (native-inputs
     `(("pkg-config" ,pkg-config)))
    (home-page "https://github.com/smarco/WFA2-lib/")
    (synopsis "Library for wavefront aligner")
    (description "The wavefront alignment (WFA) algorithm is an **exact** gap-affine algorithm that takes advantage of homologous regions between the sequences to accelerate the alignment process. Unlike to traditional dynamic programming algorithms that run in quadratic time, the WFA runs in time `O(ns+s^2)`, proportional to the sequence length `n` and the alignment score `s`, using `O(s^2)` memory (or `O(s)` using the ultralow/BiWFA mode).")
    (license license:expat))) ;; MIT license

wfa2-lib-git
