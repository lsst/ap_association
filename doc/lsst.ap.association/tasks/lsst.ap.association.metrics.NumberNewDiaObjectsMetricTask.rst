.. lsst-task-topic:: lsst.ap.association.metrics.NumberNewDiaObjectsMetricTask

#############################
NumberNewDiaObjectsMetricTask
#############################

``NumberNewDiaObjectsMetricTask`` computes the number of DIAObjects created when processing data through source association (as the ``ap_association.numNewDiaObjects`` metric).
It requires task metadata as input, and can operate at image-level or coarser granularity.

.. _lsst.ap.association.metrics.NumberNewDiaObjectsMetricTask-summary:

Processing summary
==================

``NumberNewDiaObjectsMetricTask`` reads ``AssociationTask`` statistics from task metadata associated with one or more processed images.
It uses these statistics to count the number of new DIAObjects added when processing those images.

.. _lsst.ap.association.metrics.NumberNewDiaObjectsMetricTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ap.association.metrics.NumberNewDiaObjectsMetricTask

.. _lsst.ap.association.metrics.NumberNewDiaObjectsMetricTask-butler:

Butler datasets
===============

Input datasets
--------------

:lsst-config-field:`~lsst.verify.tasks.metadataMetricTask.MetadataMetricConfig.metadata`
    The metadata of the top-level pipeline task (e.g., ``CharacterizeImageTask``, ``DiaPipeTask``) being instrumented.

.. _lsst.ap.association.metrics.NumberNewDiaObjectsMetricTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ap.association.metrics.NumberNewDiaObjectsMetricTask

.. _lsst.ap.association.metrics.NumberNewDiaObjectsMetricTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ap.association.metrics.NumberNewDiaObjectsMetricTask
