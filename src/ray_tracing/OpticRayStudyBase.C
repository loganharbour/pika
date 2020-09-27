#include "OpticRayStudyBase.h"

InputParameters
OpticRayStudyBase::validParams()
{
  auto params = RepeatableRayStudyBase::validParams();

  params.addParam<unsigned int>("energy_groups", 0, "Number of energy groups");

  params.set<bool>("_use_ray_registration") = false;
  params.set<bool>("_allow_generation_during_propagation") = true;
  params.set<bool>("_ray_dependent_subdomain_setup") = false;

  return params;
}

OpticRayStudyBase::OpticRayStudyBase(const InputParameters & parameters)
  : RepeatableRayStudyBase(parameters),
    _num_energy_groups(getParam<unsigned int>("energy_groups")),
    _energy_data_index(_num_energy_groups ? registerRayData("energy_data_g0")
                                          : Ray::INVALID_RAY_DATA_INDEX),
    _threaded_next_ids(libMesh::n_threads()),
    _threaded_rays_added(libMesh::n_threads()),
    _max_thread_ids(1000000)
{
  // We will resize the data manually when we generate the Rays
  _set_data_at_generation = false;

  if (_num_energy_groups > 1)
    for (unsigned int g = 1; g < _num_energy_groups; ++g)
      registerRayData("energy_data_g" + std::to_string(g));
}

Ray &
OpticRayStudyBase::defineRay(const Point & start,
                             const Point & direction,
                             const RayID id,
                             const Elem * elem /* = nullptr */,
                             const unsigned short incoming_side /* = invalid_side */,
                             const bool size_data /* = true */,
                             const bool size_aux_data /* = true */)
{
  std::shared_ptr<Ray> ray = _ray_pool->acquire();
  _rays.push_back(ray);

  const auto unit_direction = direction.unit();

  const auto end = rayDomainIntersection(start, unit_direction);

  ray->setStart(start);
  ray->setEnd(end);
  ray->setDirection(unit_direction);
  ray->setID(id);

  if (elem)
  {
    mooseAssert(elem->contains_point(start), "Elem doesn't contain start");
    ray->setStartingElem(elem);

    if (incoming_side != RayTracingCommon::invalid_side)
    {
      mooseAssert(elem->side_ptr(incoming_side)->contains_point(start),
                  "Side on elem doesn't contain start");
      ray->setIncomingSide(incoming_side);
    }
  }

  if (size_data)
    ray->data().resize(rayDataSize(), 0);
  if (size_aux_data)
    ray->auxData().resize(rayAuxDataSize(), 0);

  return *ray.get();
}

Ray &
OpticRayStudyBase::addReflectedRay(const Point & start,
                                   const Elem * elem,
                                   const Point & direction,
                                   const Ray & parent_ray,
                                   const Real energy_factor,
                                   const THREAD_ID tid)
{
  std::shared_ptr<Ray> ray = _ray_pool->acquire();
  addToThreadedBuffer(ray, tid);

  if (++_threaded_rays_added[tid] > _max_thread_ids)
    mooseError("bad");

  const Point end = rayDomainIntersection(start, direction);
  ray->setStartingElem(elem);
  ray->setStart(start);
  ray->setEnd(end);
  ray->setDirection(direction);
  ray->setID(_threaded_next_ids[tid]++);

  ray->auxData().resize(rayAuxDataSize(), 0);
  ray->data().resize(rayDataSize(), 0);

  for (unsigned int g = 0; g < _num_energy_groups; ++g)
    ray->data(_energy_data_index + g) = energy_factor * parent_ray.data(_energy_data_index + g);

  return *ray;
}

void
OpticRayStudyBase::buildSegmentQuadrature(const Point & start,
                                          const Point & end,
                                          std::vector<Point> & points,
                                          std::vector<Real> & weights) const
{
  // For an optic study, we only care about reinit at both ends of the segment
  // instead of the standard 1D quadrature along the segment that comes from
  // RayTracingStudy
  points = {start, end};
  weights = {1, 1};
}

void
OpticRayStudyBase::generateRays()
{
  RepeatableRayStudyBase::generateRays();

  RayID max_ray_id = 0;
  for (const auto & ray : _rays)
    max_ray_id = std::max(ray->id(), max_ray_id);
  _communicator.sum(max_ray_id);

  generateThreadedNextIDs();
}

void
OpticRayStudyBase::generateThreadedNextIDs()
{
  RayID max_ray_id = 0;
  for (const auto & ray : _rays)
    max_ray_id = std::max(ray->id(), max_ray_id);
  _communicator.max(max_ray_id);

  const RayID begin_pid =
      max_ray_id + 1 + (RayID)processor_id() * _max_thread_ids * (RayID)libMesh::n_threads();
  for (RayID tid = 0; tid < libMesh::n_threads(); ++tid)
    _threaded_next_ids[tid] = begin_pid + tid * _max_thread_ids;

  std::stringstream oss;
  oss << processor_id() << ": threaded next " << _threaded_next_ids[0] << "\n";
  std::cerr << oss.str();
}
