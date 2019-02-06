
url='guava.ufes.br/guava_repet/?data=repet&loc='
scaffold=$1
start=$2
end=$3
link=${4-$url}
pos=$scaffold%3A$start..$end

echo http://$link$pos\&tracks=genes%2Cillumina%20coverage\&highlight=$pos

