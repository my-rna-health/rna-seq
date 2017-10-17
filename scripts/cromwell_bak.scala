#! /usr/local/bin/amm

import ammonite.ops._
import ammonite.ops.ImplicitWd._

@main
def mysql(tag: String = "latest", folder: Path = pwd) = {
	%%docker("run", "-d", "-p", "3306:3306",
		"-v", s"${folder}:/var/lib/mysql",
		"-v", "./init:/docker-entrypoint-initdb.d",
		"-e", s"MYSQL_ROOT_PASSWORD=cromwell",
		"-e", s"MYSQL_DATABASE=cromwell_db",
		s"mysql:${tag}")(folder)
	//docker run --name some-mysql -v /my/own/datadir:/var/lib/mysql -e MYSQL_ROOT_PASSWORD=my-secret-pw -d mysql:tag
}

@main
def server(configName: String = "application.conf", path: Path = pwd) = {
  val p = path.toString
  //-Dconfig.file=/path/to/yourOverrides.conf cromwell.jar
  %java(s"-Dconfig.file=${p}/${configName}", "-jar", s"${p}/cromwell.jar", "server")
}
