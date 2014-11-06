#include "DiffusionLinearElement.h"
#include "utils/Utilities.h"

namespace AMP {
	namespace Operator {

		void DiffusionLinearElement::apply() {
			const std::vector<Real> & JxW = (*d_JxW);

			const std::vector<std::vector<Real> > & phi = (*d_phi);

			const std::vector<std::vector<RealGradient> > & dphi = (*d_dphi);

			std::vector<std::vector<double> > & elementStiffnessMatrix =
				(*d_elementStiffnessMatrix);

			d_fe->reinit(d_elem);

			const std::vector<Point>& q_point = d_fe->get_xyz();

			d_transportModel->preLinearElementOperation();

			const unsigned int num_local_dofs = d_elem->n_nodes();
			AMP_ASSERT(d_num_dofs == num_local_dofs);

			std::vector<double> conductivity(d_qrule->n_points());
		    std::vector<std::vector< AMP::shared_ptr<std::vector<double> > > > conductivityTensor(3,
			        std::vector<AMP::shared_ptr<std::vector<double> > >(3));
			if (d_transportModel->isaTensor()) {
				d_transportTensorModel = AMP::dynamic_pointer_cast<DiffusionTransportTensorModel>(d_transportModel);
				for (int i=0; i<3; i++) for (int j=0; j<3; j++) {
					std::vector<double> *vd = new std::vector<double>(d_qrule->n_points());
					conductivityTensor[i][j].reset(vd);
				}
			}

			if (d_transportAtGauss) {
				std::map<std::string,  AMP::shared_ptr<std::vector<double> > > args;
				if(!d_LocalTemperature.empty()) {
					std::vector<double>* temperature = new std::vector<double>(d_qrule->n_points());
					for (unsigned int qp = 0; qp < d_qrule->n_points(); qp++) {
						(*temperature)[qp] = 0.0;
						for (unsigned int j = 0; j < num_local_dofs; j++) {
							(*temperature)[qp] += d_LocalTemperature[j] * phi[j][qp];
						}//end for j
					}//end for qp
					args.insert(std::make_pair("temperature",   AMP::shared_ptr<std::vector<double> >(temperature)));
				}
				if(!d_LocalConcentration.empty()) {        
					std::vector<double>* concentration = new std::vector<double>(d_qrule->n_points());
					for (unsigned int qp = 0; qp < d_qrule->n_points(); qp++) {
						(*concentration)[qp] = 0.0;
						for (unsigned int j = 0; j < num_local_dofs; j++) {
							(*concentration)[qp] += d_LocalConcentration[j] * phi[j][qp];
						}//end for j
					}//end for qp
					args.insert(std::make_pair("concentration", AMP::shared_ptr<std::vector<double> >(concentration)));
				}
				if(!d_LocalBurnup.empty()) {
					std::vector<double>* burnup = new std::vector<double>(d_qrule->n_points());
					for (unsigned int qp = 0; qp < d_qrule->n_points(); qp++) {
						(*burnup)[qp] = 0.0;
						for (unsigned int j = 0; j < num_local_dofs; j++) {
							(*burnup)[qp] += d_LocalBurnup[j] * phi[j][qp];
						}//end for j
					}//end for qp
					args.insert(std::make_pair("burnup",        AMP::shared_ptr<std::vector<double> >(burnup)));
				}
				if (not d_transportModel->isaTensor()) {
					d_transportModel->getTransport(conductivity, args, q_point);
				} else {
					d_transportTensorModel->getTensorTransport(conductivityTensor, args, q_point);
				}

			} else {
				std::vector<double> &temperature(*new std::vector<double>(d_LocalTemperature));
				std::vector<double> &concentration(*new std::vector<double>(d_LocalConcentration));
				std::vector<double> &burnup(*new std::vector<double>(d_LocalBurnup));
				std::map<std::string,  AMP::shared_ptr<std::vector<double> > > args;
				args.insert(std::make_pair("temperature",   AMP::shared_ptr<std::vector<double> >(&temperature)));
				args.insert(std::make_pair("concentration", AMP::shared_ptr<std::vector<double> >(&concentration)));
				args.insert(std::make_pair("burnup",        AMP::shared_ptr<std::vector<double> >(&burnup)));

				std::vector<double> nodalConductivity(num_local_dofs);
			    std::vector<std::vector< AMP::shared_ptr<std::vector<double> > > > nodalConductivityTensor(3,
						std::vector<AMP::shared_ptr<std::vector<double> > >  (3));
				if (not d_transportModel->isaTensor()) {
					d_transportModel->getTransport(nodalConductivity, args, q_point);
				} else {
		        	for (int i=0; i<3; i++) for (int j=0; j<3; j++) {
		        		std::vector<double> *vec(new std::vector<double>(num_local_dofs));
		        		nodalConductivityTensor[i][j].reset(vec);
		        	}
					d_transportTensorModel->getTensorTransport(nodalConductivityTensor, args, q_point);
				}

				for (unsigned int qp = 0; qp < d_qrule->n_points(); qp++) {
					conductivity[qp] = 0.0;
					if (not d_transportModel->isaTensor()) {
						for (unsigned int n = 0; n < num_local_dofs; n++) {
							conductivity[qp] += nodalConductivity[n] * phi[n][qp];
						}//end for n
					} else {
						for (unsigned int n = 0; n < num_local_dofs; n++) {
							for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
								(*conductivityTensor[i][j])[qp] += (*nodalConductivityTensor[i][j])[n] * phi[n][qp];
							}
						}//end for n
					}
				}//end for qp
			}

			for (unsigned int qp = 0; qp < d_qrule->n_points(); qp++) {
				d_transportModel->preLinearGaussPointOperation();

				if (not d_transportModel->isaTensor()) {
					for (unsigned int n = 0; n < num_local_dofs; n++) {
						for (unsigned int k = 0; k < num_local_dofs; k++) {
							elementStiffnessMatrix[n][k] += (JxW[qp] * conductivity[qp]
									* (dphi[n][qp] * dphi[k][qp]));
						}//end for k
					}//end for n
				} else {
					for (unsigned int n = 0; n < num_local_dofs; n++) {
						for (unsigned int k = 0; k < num_local_dofs; k++) {
							for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
								elementStiffnessMatrix[n][k] += (JxW[qp] * (*conductivityTensor[i][j])[qp]
										* (dphi[n][qp](i) * dphi[k][qp](j)));
							}
						}//end for k
					}//end for n
				}

				d_transportModel->postLinearGaussPointOperation();
			}//end for qp

			d_transportModel->postLinearElementOperation();
		}

	}
}


