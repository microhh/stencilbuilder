#include <fstream>
#include <iostream>
#include <cstdio>
#include <Kokkos_Core.hpp>

namespace
{
    #ifdef SINGLE_PRECISION
    using Real = float;
    #else
    using Real = double;
    #endif
    using Real_ptr = Real* const __restrict__;

    using Array_3d_cpu = Kokkos::View<Real***, Kokkos::LayoutRight, Kokkos::HostSpace>;
    using Range_3d_cpu = Kokkos::MDRangePolicy<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::Rank<3>>;

    #ifdef KOKKOS_ENABLE_CUDA
    using Array_3d_gpu = Kokkos::View<Real***, Kokkos::LayoutRight, Kokkos::CudaSpace>;
    using Range_3d_gpu = Kokkos::MDRangePolicy<
        Kokkos::Cuda,
        Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Right>>;
    #endif

    void init(Real_ptr u, Real_ptr v, Real_ptr w, Real_ptr ut, const size_t ncells)
    {
        for (size_t i=0; i<ncells; ++i)
        {
            u [i] = pow(i,2)/pow(i+1,2);
            v [i] = pow(i,2)/pow(i+1,2);
            w [i] = pow(i,2)/pow(i+1,2);
            ut[i] = 0.;
        }
    }


    template<class Array_3d>
    struct diff
    {
        Array_3d ut;
        const Array_3d u;
        const Array_3d v;
        const Array_3d w;
        const Real visc;

        diff(
                Array_3d ut_, const Array_3d u_, const Array_3d v_, const Array_3d w_,
                const Real visc_) :
            ut(ut_), u(u_), v(v_), w(w_), visc(visc_) {}

        KOKKOS_INLINE_FUNCTION
        Real grad(const Real m2, const Real m1, const Real p1, const Real p2) const
        {
            return (1./24)*m2 + (-27./24)*m1 + (27./24.)*p1 + (-1./24.)*p2;
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(typename Array_3d::size_type k, typename Array_3d::size_type j, typename Array_3d::size_type i) const
        {
            ut(k, j, i) += visc * (
                    + grad( grad(u(k-3, j, i), u(k-2, j, i), u(k-1, j, i), u(k  , j, i)),
                            grad(u(k-2, j, i), u(k-1, j, i), u(k  , j, i), u(k+1, j, i)),
                            grad(u(k-1, j, i), u(k  , j, i), u(k+1, j, i), u(k+2, j, i)),
                            grad(u(k  , j, i), u(k+1, j, i), u(k+2, j, i), u(k+3, j, i)) )

                    + grad( grad(u(k, j-3, i), u(k, j-2, i), u(k, j-1, i), u(k, j  , i)),
                            grad(u(k, j-2, i), u(k, j-1, i), u(k, j  , i), u(k, j+1, i)),
                            grad(u(k, j-1, i), u(k, j  , i), u(k, j+1, i), u(k, j+2, i)),
                            grad(u(k, j  , i), u(k, j+1, i), u(k, j+2, i), u(k, j+3, i)) )

                    + grad( grad(u(k, j, i-3), u(k, j, i-2), u(k, j, i-1), u(k, j, i  )),
                            grad(u(k, j, i-2), u(k, j, i-1), u(k, j, i  ), u(k, j, i+1)),
                            grad(u(k, j, i-1), u(k, j, i  ), u(k, j, i+1), u(k, j, i+2)),
                            grad(u(k, j, i  ), u(k, j, i+1), u(k, j, i+2), u(k, j, i+3)) ) );
        }
    };
}

int main(int argc, char* argv[])
{
    Kokkos::initialize(argc, argv);
    {
        if (argc != 2)
        {
            std::cout << "Add the grid size as an argument!" << std::endl;
            return 1;
        }

        constexpr int nloop = 30;
        const size_t itot = std::stoi(argv[1]);
        const size_t jtot = std::stoi(argv[1]);
        const size_t ktot = std::stoi(argv[1]);
        const size_t ncells = itot*jtot*ktot;

        constexpr Real visc = 0.1;
        constexpr Real dxidxi = 0.1;
        constexpr Real dyidyi = 0.1;
        constexpr Real dzidzi = 0.1;

        Kokkos::Timer timer;

        #ifdef KOKKOS_ENABLE_CUDA
        // SOLVE ON THE GPU.
        Array_3d_gpu a_gpu ("a_gpu" , ktot, jtot, itot);
        Array_3d_gpu at_gpu("at_gpu", ktot, jtot, itot);

        Range_3d_gpu range_3d_gpu({1, 1, 1}, {ktot-1, jtot-1, itot-1}, {1, 1, 64});

        Array_3d_gpu::HostMirror a_tmp = Kokkos::create_mirror_view(a_gpu);
        Array_3d_gpu::HostMirror at_tmp = Kokkos::create_mirror_view(at_gpu);

        init(a_tmp.data(), at_tmp.data(), ncells);

        Kokkos::deep_copy(a_gpu, a_tmp);
        Kokkos::deep_copy(at_gpu, at_tmp);

        // Time performance.
        timer.reset();

        for (int i=0; i<nloop; ++i)
        {
            Kokkos::parallel_for(
                    range_3d_gpu,
                    diff<Array_3d_gpu>(at_gpu, a_gpu, visc, dxidxi, dyidyi, dzidzi));
        }

        Kokkos::fence();

        double duration_gpu = timer.seconds();

        printf("time/iter (GPU) = %E s (%i iters)\n", duration_gpu/(double)nloop, nloop);

        Kokkos::deep_copy(at_tmp, at_gpu);

        printf("at=%.20f\n", at_tmp.data()[itot*jtot+itot+itot/4]);
        #endif

        // SOLVE ON THE CPU.
        Array_3d_cpu ut_cpu("ut_cpu", ktot+8, jtot+8, itot+8);
        Array_3d_cpu u_cpu ("u_cpu" , ktot+8, jtot+8, itot+8);
        Array_3d_cpu v_cpu ("v_cpu" , ktot+8, jtot+8, itot+8);
        Array_3d_cpu w_cpu ("w_cpu" , ktot+8, jtot+8, itot+8);

        Range_3d_cpu range_3d_cpu({4, 4, 4}, {ktot-4, jtot-4, itot-4}, {2, 0, 0});

        init(u_cpu.data(), v_cpu.data(), w_cpu.data(), ut_cpu.data(), ncells);

        // Time performance.
        timer.reset();

        for (int i=0; i<nloop; ++i)
        {
            Kokkos::parallel_for(
                    range_3d_cpu,
                    diff<Array_3d_cpu>(ut_cpu, u_cpu, v_cpu, w_cpu, visc));
        }

        Kokkos::fence();

        double duration_cpu = timer.seconds();

        printf("time/iter (CPU) = %E s (%i iters)\n", duration_cpu/(double)nloop, nloop);

        printf("at=%.20f\n", ut_cpu.data()[itot*jtot+itot+itot/4]);

        /*
        std::ofstream binary_file("at_kokkos_cpu.bin", std::ios::out | std::ios::trunc | std::ios::binary);

        if (binary_file)
            binary_file.write(reinterpret_cast<const char*>(at_cpu.data()), ncells*sizeof(Real));
        else
        {
            std::string error = "Cannot write file \"at_cuda.bin\"";
            throw std::runtime_error(error);
        }
        */
    }

    Kokkos::finalize();

    return 0;
}
