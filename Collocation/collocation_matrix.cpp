#include <container.h>
#include <collocation_matrix.h>

void Container::subtract_delta(const Collocation_Matrix &colm, double coeff){
	for(int i = 0; i < colm.sz; ++i)
		data[i] -= colm.delta(i);
	return;
}