# golob/count_combined:0.2
FROM golob/fastatools:0.4__bcw.0.3.0

ADD count_combined.py /usr/local/bin/
RUN chmod +x /usr/local/bin/count_combined.py
RUN apk add --no-cache python3-dev build-base
RUN pip3 install pandas
RUN pip3 install biopython

