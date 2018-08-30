from __future__ import absolute_import, division
from __future__ import print_function
import logging
logger = logging.getLogger(__name__)

from dials.util import log

debug_handle = log.debug_handle(logger)
info_handle = log.info_handle(logger)

from libtbx.utils import Sorry

from dials_algorithms_indexing_ext import *

import iotbx.phil # implicit import

from dials.array_family import flex
from cctbx import crystal
from dials.algorithms.indexing.stills_indexer import stills_indexer



class iota_refiner(stills_indexer):
  ''' Class for doing refinement and outlier rejection after iota indexing is done
      Although this class subclasses stills_indexer, it should never be used for any
      iota related indexing. There is a separate class for that.'''
  def __init__ (self, experiments, indexed, imagesets,params):
    self.all_params = params
    self.reflections = indexed
    self.experiments = experiments
    self.imagesets = imagesets

  def run_refinement_and_outlier_rejection(self):
    ''' Code taken from the index function of stills_indexer '''
    self.d_min = self.all_params.indexing.refinement_protocol.d_min_start
    self.indexed_reflections = (self.reflections['id'] > -1)
    if self.d_min is None:
      sel = self.reflections['id'] <= -1
    else:
      sel = flex.bool(len(self.reflections), False)
      lengths = 1/self.reflections['rlp'].norms()
      isel = (lengths >= self.d_min).iselection()
      sel.set_selected(isel, True)
      sel.set_selected(self.reflections['id'] > -1, False)
    self.unindexed_reflections = self.reflections.select(sel)

    reflections_for_refinement = self.reflections.select(
      self.indexed_reflections)
    import copy
    experiments = copy.deepcopy(self.experiments)
    print("Starting Refinement")

    try:
      refined_experiments, refined_reflections = self.refine(
        experiments, reflections_for_refinement)
    except Exception as e:
      s = str(e)
      raise Sorry(e)

    # sanity check for unrealistic unit cell volume increase during refinement
    # usually this indicates too many parameters are being refined given the
    # number of observations provided.
    if not self.all_params.indexing.refinement_protocol.disable_unit_cell_volume_sanity_check:
      for orig_expt, refined_expt in zip(experiments, refined_experiments):
        uc1 = orig_expt.crystal.get_unit_cell()
        uc2 = refined_expt.crystal.get_unit_cell()
        volume_change = abs(uc1.volume()-uc2.volume())/uc1.volume()
        cutoff = 0.5
        if volume_change > cutoff:
          msg = "\n".join((
            "Unrealistic unit cell volume increase during refinement of %.1f%%.",
            "Please try refining fewer parameters, either by enforcing symmetry",
            "constraints (space_group=) and/or disabling experimental geometry",
            "refinement (detector.fix=all and beam.fix=all). To disable this",
            "sanity check set disable_unit_cell_volume_sanity_check=True.")) %(
            100*volume_change)
          raise Sorry(msg)

    self.refined_reflections = refined_reflections.select(
      refined_reflections['id'] > -1)

    for i, imageset in enumerate(self.imagesets):
      ref_sel = self.refined_reflections.select(
        self.refined_reflections['imageset_id'] == i)
      ref_sel = ref_sel.select(ref_sel['id'] >= 0)
      for i_expt in set(ref_sel['id']):
        expt = refined_experiments[i_expt]
        imageset.set_detector(expt.detector)
        imageset.set_beam(expt.beam)
        imageset.set_goniometer(expt.goniometer)
        imageset.set_scan(expt.scan)
        expt.imageset = imageset


    if not (self.all_params.refinement.parameterisation.beam.fix == 'all'
            and self.all_params.refinement.parameterisation.detector.fix == 'all'):
      # Experimental geometry may have changed - re-map centroids to
      # reciprocal space

      spots_mm = self.reflections
      self.reflections = flex.reflection_table()
      for i, imageset in enumerate(self.imagesets):
        spots_sel = spots_mm.select(spots_mm['imageset_id'] == i)
        self.map_centroids_to_reciprocal_space(
          spots_sel, imageset.get_detector(), imageset.get_beam(),
          imageset.get_goniometer())
        self.reflections.extend(spots_sel)

    # update for next cycle
    experiments = refined_experiments
    self.refined_experiments = refined_experiments

    # discard experiments with zero reflections after refinement
    id_set = set(self.refined_reflections['id'])
    if len(id_set) < len(self.refined_experiments):
      filtered_refined_reflections = flex.reflection_table()
      for i in xrange(len(self.refined_experiments)):
        if i not in id_set:
          del self.refined_experiments[i]
      for old, new in zip(sorted(id_set), range(len(id_set))):
        subset = self.refined_reflections.select(self.refined_reflections['id'] == old)
        subset['id'] = flex.int(len(subset), new)
        filtered_refined_reflections.extend(subset)
      self.refined_reflections = filtered_refined_reflections

    if len(self.refined_experiments) > 1:
      from dials.algorithms.indexing.compare_orientation_matrices \
           import show_rotation_matrix_differences
      # FIXME
      #show_rotation_matrix_differences(
      #  self.refined_experiments.crystals(), out=info_handle)

    #logger.info("Final refined crystal models:")
    for i, crystal_model in enumerate(self.refined_experiments.crystals()):
      n_indexed = 0
      for i_expt in experiments.where(crystal=crystal_model):
        n_indexed += (self.reflections['id'] == i).count(True)
      #logger.info("model %i (%i reflections):" %(i+1, n_indexed))
      #logger.info(crystal_model)

    if 'xyzcal.mm' in self.refined_reflections: # won't be there if refine_all_candidates = False and no isoforms
      self.refined_reflections['xyzcal.px'] = flex.vec3_double(
        len(self.refined_reflections))
      for i, imageset in enumerate(self.imagesets):
        imgset_sel = self.refined_reflections['imageset_id'] == i
        # set xyzcal.px field in self.refined_reflections
        refined_reflections = self.refined_reflections.select(imgset_sel)
        panel_numbers = flex.size_t(refined_reflections['panel'])
        xyzcal_mm = refined_reflections['xyzcal.mm']
        x_mm, y_mm, z_rad = xyzcal_mm.parts()
        xy_cal_mm = flex.vec2_double(x_mm, y_mm)
        xy_cal_px = flex.vec2_double(len(xy_cal_mm))
        for i_panel in range(len(imageset.get_detector())):
          panel = imageset.get_detector()[i_panel]
          sel = (panel_numbers == i_panel)
          isel = sel.iselection()
          ref_panel = refined_reflections.select(panel_numbers == i_panel)
          xy_cal_px.set_selected(
            sel, panel.millimeter_to_pixel(xy_cal_mm.select(sel)))
        x_px, y_px = xy_cal_px.parts()
        scan = imageset.get_scan()
        if scan is not None:
          z_px = scan.get_array_index_from_angle(z_rad, deg=False)
        else:
          # must be a still image, z centroid not meaningful
          z_px = z_rad
        xyzcal_px = flex.vec3_double(x_px, y_px, z_px)
        self.refined_reflections['xyzcal.px'].set_selected(imgset_sel, xyzcal_px)

        return self.refined_experiments, self.refined_reflections
