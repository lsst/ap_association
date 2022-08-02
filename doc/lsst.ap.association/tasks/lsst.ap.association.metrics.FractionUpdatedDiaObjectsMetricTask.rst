.. lsst-task-topic:: lsst.ap.association.metrics.FractionUpdatedDiaObjectsMetricTask

###################################
FractionUpdatedDiaObjectsMetricTask
###################################

``FractionUpdatedDiaObjectsMetricTask`` computes the fraction of DIAObjects that were updated when processing data through source association (as the ``ap_association.fracUpdatedDiaObjects`` metric).
It requires task metadata as input.
While the task can operate at image-level or coarser granularity, the current algorithm may double-count objects and should not be run on multiple visits.

.. _lsst.ap.association.metrics.FractionUpdatedDiaObjectsMetricTask-summary:

Processing summary
==================

``FractionUpdatedDiaObjectsMetricTask`` reads ``AssociationTask`` statistics from task metadata associated with one or more processed images.
It uses these statistics to compute the fraction of potentially updatable DIAObjects that were updated with new sources when processing those images.

.. _lsst.ap.association.metrics.FractionUpdatedDiaObjectsMetricTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ap.association.metrics.FractionUpdatedDiaObjectsMetricTask

.. _lsst.ap.association.metrics.FractionUpdatedDiaObjectsMetricTask-butler:

Butler datasets
===============

Input datasets
--------------

:lsst-config-field:`~lsst.verify.tasks.metadataMetricTask.MetadataMetricConfig.metadata`
    The metadata of the top-level pipeline task (e.g., ``CharacterizeImageTask``, ``DiaPipeTask``) being instrumented.

.. _lsst.ap.association.metrics.FractionUpdatedDiaObjectsMetricTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ap.association.metrics.FractionUpdatedDiaObjectsMetricTask

.. _lsst.ap.association.metrics.FractionUpdatedDiaObjectsMetricTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ap.association.metrics.FractionUpdatedDiaObjectsMetricTask
