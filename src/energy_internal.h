/*
 * energy_internal.h —— 仅供测试使用的内部接口，不属于公共 API。
 *
 * 负荷 L 与活动量 A 由库根据心率 / 步数 / 热量计算，调用方不传。
 * 单元测试需要在固定的 L / A 下核对动力学与观测方程，因此保留下面两个带覆盖值的入口：
 * load / activity 为有限且 >= 0 时直接采用（裁剪到 [0,1]），否则由库计算。
 */
#ifndef ENERGY_INTERNAL_H
#define ENERGY_INTERNAL_H

#include "energy.h"

#ifdef __cplusplus
extern "C" {
#endif

int energy_ctx_step_fixed(energy_ctx_t *ctx, const energy_sample_t *sample,
                          double load, double activity, energy_point_t *out);

/* fix_load / fix_act 为长度 n 的数组或 NULL；单个元素 < 0 表示该点由库计算。 */
int energy_compute_fixed(const energy_sample_t *samples, size_t n,
                         const energy_config_t *cfg,
                         const double *fix_load, const double *fix_act,
                         energy_output_t *out);

#ifdef __cplusplus
}
#endif

#endif /* ENERGY_INTERNAL_H */
