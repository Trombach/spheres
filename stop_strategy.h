#ifndef STOP_STRATEGY
#define STOP_STRATEGY

#include <iostream>
#include <iomanip>
#include <dlib/matrix.h>
#include <dlib/algs.h>

namespace dlib {

	class stop_strategy {
		public:
			explicit stop_strategy (
				double min_delta = 1e-7
			) : _verbose(false), _been_used(false), _min_delta(min_delta), _max_iter(0), _cur_iter(0), _prev_funct_value(0) 
			{
				DLIB_ASSERT (
					min_delta >= 0,
					"\t stop_strategy(min_delta)"
					<< "\n\t min_delta can't be negative"
					<< "\n\t min_delta: " << min_delta
					);
			}

			stop_strategy (
				double min_delta,
				unsigned long max_iter
			) : _verbose(false), _been_used(false), _min_delta(min_delta), _max_iter(max_iter), _cur_iter(0), _prev_funct_value(0)
			{
				DLIB_ASSERT (
					min_delta >= 0,
					"\t stop_strategy(min_delta)"
					<< "\n\t min_delta can't be negative"
					<< "\n\t min_delta: " << min_delta
					<< "\n\t max_iter: " << max_iter
					);
			}

			stop_strategy& be_verbose(std::ostream &out)
			{
				_verbose = true;
				output = &out;
				output->precision(16);
				*output << std::scientific;
				return *this;
			}

			template <typename T>
			bool should_continue_search (
				const T& ,
				const double funct_value,
				const T& funct_derivative
			)
			{
				if (_verbose)
				{
					*output << std::setw(4) << std::setfill(' ') << std::left;
					*output << _cur_iter 
						<< std::setw(25) << funct_value << length(funct_derivative) << std::endl;
				}

				++_cur_iter;
				if (_been_used)
				{
					if (_max_iter != 0 && _cur_iter > _max_iter)
					{
						*output << "nsteps reached" << std::endl;
						return false;
					}

					if (std::abs(funct_value - _prev_funct_value) < _min_delta)
					{
						*output << "Stationary point found" << std::endl;
						return false;
					}
				}

				_been_used = true;
				_prev_funct_value = funct_value;
				return true;
			}


		private:
			bool _verbose;

			bool _been_used;
			double _min_delta;
			unsigned long _max_iter;
			unsigned long _cur_iter;
			double _prev_funct_value;

			std::ostream *output;
	};
}

#endif
