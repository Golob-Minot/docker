## VERISON v1.1.alpha19-4-g1189285

FROM ubuntu:trusty
RUN mkdir /fh && mkdir /app
RUN mkdir /pplacer && mkdir /pplacer/src && mkdir /data
WORKDIR /pplacer/src
RUN apt-get -y update && apt-get -y install camlp4-extra gawk libgsl0-dev libsqlite3-dev libz-dev m4 wget make ocaml patch build-essential git pkg-config unzip python-pip python-dev sqlite3
RUN wget https://raw.github.com/ocaml/opam/master/shell/opam_installer.sh -O - | sh -s /usr/local/bin
RUN ln -s /usr/local/bin/opam /usr/bin/opam && /usr/local/bin/opam init -y
#RUN command -v ocamlc && ( ocamlc -version | grep -q 3.12.1 ) || opam switch install 3.12.1
RUN opam repo add pplacer-deps http://matsen.github.com/pplacer-opam-repository &&  opam update pplacer-deps && eval `opam config env`
RUN git clone https://github.com/matsen/pplacer
WORKDIR /pplacer/src/pplacer
#RUN opam install -y depext && opam depext -y csv.1.5 ounit.2.0.0 xmlm.1.1.1 mcl.12-068oasis4 batteries.2.7.0 ocaml-gsl.0.6.3 sqlite3.4.1.3 camlzip.1.05 && opam install -y csv.1.5 ounit.2.0.0 xmlm.1.1.1 mcl.12-068oasis4 batteries.2.7.0 ocaml-gsl.0.6.3 sqlite3.4.1.3 camlzip.1.05
RUN opam install -y depext && opam depext -y csv ounit xmlm mcl batteries ocaml-gsl sqlite3 camlzip && opam install -y csv ounit xmlm mcl batteries ocaml-gsl sqlite3 camlzip
RUN eval $(opam config env) && make
RUN cp /pplacer/src/pplacer/bin/* /usr/local/bin
RUN pip install pip --upgrade && pip install numpy biopython pandas --upgrade
WORKDIR /pplacer/src/pplacer/scripts
RUN python setup.py install
RUN pip install taxtastic==0.8.3
#cleanup
WORKDIR /pplacer
RUN pip install csvkit==1.0.2
