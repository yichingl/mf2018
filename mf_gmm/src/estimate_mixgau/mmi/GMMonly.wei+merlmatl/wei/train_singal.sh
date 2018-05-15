cp ~/projects/dtrain/optimized/src/gmmtr bin/gmmtr

for i in $(seq 4 4); do
./bin/gmmtr -i proto/g.d2.m$i.f.pro -l script/g.n1.d2.m4.f.tr1000 -o model/g.n1.d2.m$i.f.tr1000.def -i1 10 -i2 20
done
