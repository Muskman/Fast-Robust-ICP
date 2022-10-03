#include <iostream>
#include "ICP.h"
#include "io_pc.h"
#include "FRICP.h"

int main(int argc, char const ** argv)
{
    typedef double Scalar;
    typedef Eigen::Matrix<Scalar, 3, Eigen::Dynamic> Vertices;
    typedef Eigen::Matrix<Scalar, 3, 1> VectorN;
    std::string folder_data;
    std::string file_source;
    std::string file_target;
    std::string file_init;
    std::string file_gt;
    std::string file_param;
    std::string res_trans_path;
    std::string out_path;
    bool use_init = false;
    bool has_groundtruth = true;
    bool print_output = true;
    MatrixXX res_trans;
    MatrixXX gt_trans;
    enum Method{ICP, AA_ICP, FICP, RICP, SparseICP, StochasticICP, SpiderICP, PPL, RPPL, SICPPPL} method=RICP;
    if(argc == 3)
    {
        folder_data = argv[1];
        file_target = folder_data + "/target.ply";
        file_source = folder_data + "/source.ply";
        file_init = folder_data + "/init.txt";
        file_gt = folder_data + "/gt.txt";
        file_param = folder_data + "/params.txt";
        out_path = folder_data + "/res/";
        method = Method(std::stoi(argv[2]));
    }
    else if(argc==2)
    {
        file_target = argv[1];
        file_source = argv[2];
        out_path = argv[3];
    }
    else
    {
        std::cout << "Usage: target.ply source.ply out_path <Method>" << std::endl;
        std::cout << "Method :\n"
                  << "0: ICP\n1: AA-ICP\n2: Fast ICP\n3: Robust ICP\n4: ICP Point-to-plane\n"
                  << "5: Robust ICP point to plane\n6: Sparse ICP\n7: Sparse ICP point to plane\n" 
                  << "8: Stochastic ICP\n9: Spider ICP" << std::endl;
        exit(0);
    }
    int dim = 3;

    for (int m=0; m<7; m++) {
        method = Method(m);

        //--- Model that will be rigidly transformed
        Vertices vertices_source, normal_source, src_vert_colors;
        read_file(vertices_source, normal_source, src_vert_colors, file_source);
        std::cout << "source: " << vertices_source.rows() << "x" << vertices_source.cols() << std::endl;

        //--- Model that source will be aligned to
        Vertices vertices_target, normal_target, tar_vert_colors;
        read_file(vertices_target, normal_target, tar_vert_colors, file_target);
        std::cout << "target: " << vertices_target.rows() << "x" << vertices_target.cols() << std::endl;

        // scaling
        Eigen::Vector3d source_scale, target_scale;
        source_scale = vertices_source.rowwise().maxCoeff() - vertices_source.rowwise().minCoeff();
        target_scale = vertices_target.rowwise().maxCoeff() - vertices_target.rowwise().minCoeff();
        double scale = std::max(source_scale.norm(), target_scale.norm());
        std::cout << "scale = " << scale << std::endl;
        vertices_source /= scale;
        vertices_target /= scale;

        /// De-mean
        VectorN source_mean, target_mean;
        source_mean = vertices_source.rowwise().sum() / double(vertices_source.cols());
        target_mean = vertices_target.rowwise().sum() / double(vertices_target.cols());
        vertices_source.colwise() -= source_mean;
        vertices_target.colwise() -= target_mean;

        double time;
        // set ICP parameters
        ICP::Parameters pars;

        // set Sparse-ICP parameters
        SICP::Parameters spars;
        spars.p = 0.4;
        spars.print_icpn = false;

        /// Initial transformation
        if(use_init)
        {
            MatrixXX init_trans;
            read_transMat(init_trans, file_init);
            init_trans.block(0, dim, dim, 1) /= scale;
            init_trans.block(0,3,3,1) += init_trans.block(0,0,3,3)*source_mean - target_mean;
            pars.use_init = true;
            pars.init_trans = init_trans;
            spars.init_trans = init_trans;
        }

        /// GT transformation
        if(has_groundtruth)
        {
            MatrixXX gt_trans;
            read_transMat(gt_trans, file_gt);
            gt_trans.block(0, dim, dim, 1) /= scale;
            gt_trans.block(0,3,3,1) += gt_trans.block(0,0,3,3)*source_mean - target_mean;
            pars.has_groundtruth = true;
            pars.gt_trans = gt_trans;
            spars.gt_trans = gt_trans;
        }

        /// print output to file
        if(print_output)
        {
            pars.print_output = print_output;
            pars.out_path = out_path;
            spars.print_output = print_output;
            spars.out_path = out_path;
        }


        ///--- Execute registration
        std::cout << "begin registration..." << std::endl;
        FRICP<3> fricp;
        double begin_reg = omp_get_wtime();
        double converge_rmse = 0;
        switch(method)
        {
        case ICP:
        {
            std::cout << "ICP" << std::endl;
            pars.f = ICP::NONE;
            pars.use_AA = false;
            pars.out_path = out_path + "m" + std::to_string(method) + "output.txt";
            fricp.point_to_point(vertices_source, vertices_target, source_mean, target_mean, pars);
            res_trans = pars.res_trans;
            break;
        }
        case AA_ICP:
        {
            std::cout << "AA-ICP" << std::endl;
            pars.out_path = out_path + "m" + std::to_string(method) + "output.txt";
            AAICP::point_to_point_aaicp(vertices_source, vertices_target, source_mean, target_mean, pars);
            res_trans = pars.res_trans;
            break;
        }
        case FICP:
        {
            std::cout << "FICP" << std::endl;
            pars.f = ICP::NONE;
            pars.use_AA = true;
            pars.out_path = out_path + "m" + std::to_string(method) + "output.txt";
            fricp.point_to_point(vertices_source, vertices_target, source_mean, target_mean, pars);
            res_trans = pars.res_trans;
            break;
        }
        case RICP:
        {
            std::cout << "FRICP" << std::endl;
            pars.f = ICP::WELSCH;
            pars.use_AA = true;
            pars.out_path = out_path + "m" + std::to_string(method) + "output.txt";
            fricp.point_to_point(vertices_source, vertices_target, source_mean, target_mean, pars);
            res_trans = pars.res_trans;
            break;
        }
        case PPL:
        {
            std::cout << "ICP PPL" << std::endl;
            pars.f = ICP::NONE;
            pars.use_AA = false;
            pars.out_path = out_path + "m" + std::to_string(method) + "output.txt";
            if(normal_target.size()==0)
            {
                std::cout << "Warning! The target model without normals can't run Point-to-plane method!" << std::endl;
                exit(0);
            }
            fricp.point_to_plane(vertices_source, vertices_target, normal_source, normal_target, source_mean, target_mean, pars);
            res_trans = pars.res_trans;
            break;
        }
        case RPPL:
        {
            std::cout << "FRICP PPL" << std::endl;
            pars.nu_end_k = 1.0/6;
            pars.f = ICP::WELSCH;
            pars.use_AA = true;
            pars.out_path = out_path + "m" + std::to_string(method) + "output.txt";
            if(normal_target.size()==0)
            {
                std::cout << "Warning! The target model without normals can't run Point-to-plane method!" << std::endl;
                exit(0);
            }
            fricp.point_to_plane_GN(vertices_source, vertices_target, normal_source, normal_target, source_mean, target_mean, pars);
            res_trans = pars.res_trans;
            break;
        }
        case SparseICP:
        {
            std::cout << "Sparse ICP" << std::endl;
            pars.out_path = out_path + "m" + std::to_string(method) + "output.txt";
            SICP::point_to_point(vertices_source, vertices_target, source_mean, target_mean, spars);
            res_trans = spars.res_trans;
            break;
        }
        case SICPPPL:
        {
            std::cout << "Sparse ICP PPL" << std::endl;
            pars.out_path = out_path + "m" + std::to_string(method) + "output.txt";
            if(normal_target.size()==0)
            {
                std::cout << "Warning! The target model without normals can't run Point-to-plane method!" << std::endl;
                exit(0);
            }
            SICP::point_to_plane(vertices_source, vertices_target, normal_target, source_mean, target_mean, spars);
            res_trans = spars.res_trans;
            break;
        }
        case StochasticICP:
        {
            std::cout << "Stochastic ICP" << std::endl;
            pars.f = ICP::NONE;
            pars.max_icp = 1000;
            pars.use_AA = false;
            pars.use_stochastic = true;
            pars.batch_ratio = 0.1;
            pars.out_path = out_path + "m" + std::to_string(method) + "output.txt";
            fricp.point_to_point_stochastic(vertices_source, vertices_target, source_mean, target_mean, pars);
            res_trans = pars.res_trans;
            break;
        }
        case SpiderICP:
        {
            std::cout << "Spider ICP" << std::endl;
            pars.f = ICP::NONE;
            pars.max_icp = 1000;
            pars.use_AA = false;
            pars.use_stochastic = true;
            pars.batch_ratio = 0.1;
            pars.inner_batch_ratio = 0.01;
            pars.q_spider = 3;
            pars.out_path = out_path + "m" + std::to_string(method) + "output.txt";
            fricp.point_to_point_spider(vertices_source, vertices_target, source_mean, target_mean, pars);
            res_trans = pars.res_trans;
            break;
        }
        }
        std::cout << "Registration done!" << std::endl;
        double end_reg = omp_get_wtime();
        time = end_reg - begin_reg;
        vertices_source = scale * vertices_source;

        Eigen::Affine3d res_T;
        res_T.linear() = res_trans.block(0,0,3,3);
        res_T.translation() = res_trans.block(0,3,3,1);
        res_trans_path = out_path + "m" + std::to_string(method) + "trans.txt";
        std::ofstream out_trans(res_trans_path);
        res_trans.block(0,3,3,1) *= scale;
        out_trans << res_trans << std::endl;
        out_trans.close();

        ///--- Write result to file
        std::string file_source_reg = out_path + "m" + std::to_string(method) + "reg_pc.ply";
        write_file(file_source, vertices_source, normal_source, src_vert_colors, file_source_reg);

        std::cout << "-----------------------------------" << std::endl;
    }

    return 0;
}
