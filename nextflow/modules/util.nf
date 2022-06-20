import java.nio.file.Paths

def verify_filepath ( file_prefixes, file_suffixes, file_extensions ) {

    def file_path = []

    for (suff in file_suffixes){
        for(ext in file_extensions){
          if ( file_prefixes ) {
            for (prefix in file_prefixes) {
              
              dir_glob = params.directory.replaceAll(/\/+$/, "") + '**'
              file_glob = prefix + suff + ext
              search_path = Paths.get(dir_glob, file_glob )
              file_path.add(search_path.toString())

              }
          } else {

              dir_glob = params.directory.replaceAll(/\/+$/, "") + '**'
              file_glob = suff + ext
              search_path = Paths.get(dir_glob, file_glob)
              file_path.add(search_path.toString())
          }
        }
    }

    return file_path
}

