FROM mysql:5.7

COPY init_user.sql /docker-entrypoint-initdb.d

ENV MYSQL_ROOT_PASSWORD=cromwell
ENV MYSQL_DATABASE=cromwell_db

EXPOSE 3306
CMD ["mysqld"]